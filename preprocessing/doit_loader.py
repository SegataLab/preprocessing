#!/usr/bin/env python
#Author: Duy Tin Truong (duytin.truong@unitn.it)
#        at CIBIO, University of Trento, Italy

__author__  = 'Duy Tin Truong (duytin.truong@unitn.it)'
__version__ = '0.1'
__date__    = '11 June 2015'

import sys
import os
import argparse 
import subprocess
import copy

from doit.task import dict_to_task, Task
from doit.cmd_base import TaskLoader
from doit.doit_cmd import DoitMain
from doit.task import clean_targets


class DoitLoader(TaskLoader):
    task_list = []
    config = {'verbosity':2}
    task_id_counter = 0

    @staticmethod
    def load_tasks(action, opt_values, pos_args):
        return DoitLoader.task_list, DoitLoader.config

    @staticmethod
    def add_task(
                 targets, 
                 file_dep, 
                 actions, 
                 name=None, 
                 pipe=False,
                 **kwargs):
        if type(targets) not in [list, tuple]:
            targets = [targets]
        if type(actions) not in [list, tuple]:
            actions = [actions]
        if type(file_dep) not in [list, tuple]:
            file_dep = [file_dep]
        if name == None:
            DoitLoader.task_id_counter += 1
            name = 'task_%d'%(DoitLoader.task_id_counter)
        actions = DoitLoader.format_actions(actions)
        if pipe:
            for i in range(len(actions)-1):
                actions[i][1]['stdout'] = subprocess.PIPE
        params = [{'name':'actions', 'default':copy.deepcopy(actions)},
                  {'name':'pipe', 'default':pipe}]

        def run(actions):
            for i in range(len(actions)):
                for k in ['stdin', 'stdout', 'stderr']:
                    if k in actions[i][1] and isinstance(actions[i][1][k], str):
                        mode = 'r' if k == 'stdin' else 'w'
                        actions[i][1][k] = open(actions[i][1][k], mode)

            finished = True
            for i in range(len(actions)):
                if actions[i][0].__class__.__name__ == 'function':
                    returncode = actions[i][0](**actions[i][1])
                    finished = finished and returncode
                else:
                    p = subprocess.Popen(
                                         actions[i][0],
                                         **actions[i][1]
                                        )
                    connect = False
                    if 'stdout' in actions[i][1]:
                        if actions[i][1]['stdout'] == subprocess.PIPE \
                            and i < len(actions) - 1:
                            connect = True
                    if connect:
                        actions[i+1][1]['stdin'] = p.stdout
                    else:
                        p.communicate()
                        finished = finished and p.returncode != None

            for i in range(len(actions)):
                for k in ['stdin', 'stdout', 'stderr']:
                    if k in actions[i][1] and isinstance(actions[i][1][k], file):
                        actions[i][1][k].close() 
            return finished


        task = {
                'name': name,
                'actions': [(run, [actions])],
                'targets': targets,
                'file_dep': file_dep,
                'params': params,
                'clean': [clean_targets],
                'title': DoitLoader.print_action                
               }
        for k, v in kwargs.items():
            task[k] = v
        DoitLoader.task_list.append(dict_to_task(task))

    @staticmethod
    def add_task2(task):
        DoitLoader.task_list.append(dict_to_task(task))

    @staticmethod
    def format_actions(actions):
        if len(actions) == 2 and \
            actions[0].__class__.__name__ in ['list', 'tuple', 'str', 'function'] and \
            actions[1].__class__.__name__ == 'dict':
            actions = [actions]
        actions = [action for action in actions if 
                    (action.__class__.__name__ == 'function') or len(action)]
        for i in range(len(actions)):
            if isinstance(actions[i], str): 
                actions[i] = [actions[i].split(), {}]
            elif actions[i].__class__.__name__ == 'function':
                actions[i] = [actions[i], {}]
            elif isinstance(actions[i], (list, tuple)):
                if len(actions[i]) == 1 or not isinstance(actions[i][1], dict):
                    if isinstance(actions[i][0], str):
                        actions[i] = [actions[i][0].split(), {}]
                    else:
                        actions[i] = [actions[i], {}]
                elif isinstance(actions[i][0], str):
                    actions[i][0] = actions[i][0].split()
        return actions

    @staticmethod
    def params2struct(params):
        result = {} 
        for p in params:
            result[p['name']] = p['default']
        result = type('Struct', (object,), result)()
        return result

    @staticmethod
    def print_action(task):
        action = 'Task %s: '%task.name
        if task.actions[0].__class__.__name__ == 'PythonAction':
            params = DoitLoader.params2struct(task.params)
            actions = params.actions   
            for i in range(len(actions)):
                connect = False
                if actions[i][0].__class__.__name__ == 'function':
                    action += str(actions[i][0])
                else:
                    if 'stdin' in actions[i][1]:
                        action += '%s > '%actions[i][1]['stdin']
                    action += ' '.join(actions[i][0])
                    if 'stdout' in actions[i][1]:
                        if actions[i][1]['stdout'] == subprocess.PIPE:
                            connect = True
                        else:
                            action += ' > %s'%actions[i][1]['stdout']
                if len(actions[i][1]):
                    action += ', ' + str(actions[i][1])

                if connect:
                    action += ' | '
                else:
                    action += ' ; '
            action = action[:-3]
        else:
            action += str(task.actions)
        return action
    
    @staticmethod
    def run(doit_args):
        code = DoitMain(DoitLoader()).run(doit_args)
        if code == 0:
            print 'Doit: all tasks were executed.'
        elif code == 1:
            print 'Doit: some tasks were failed.'
        elif code == 2:
            print 'Doit: error when executing tasks.'
        elif code == 3:
            print 'Doit: error before task execution starts.'
        return code


def read_params():
    p = argparse.ArgumentParser()
    p.add_argument(
        '--doit_db', 
        required=False, 
        default='.doit.db',
        type=str,
        help='Doit database.')
    p.add_argument(
        '--nprocs', 
        required=False, 
        default=1, 
        type=int,
        help='Number of processors.')
    p.add_argument(
        '--use_threads', 
        required=False, 
        dest='use_threads',
        action='store_true',
        help='Use threads.')
    p.set_defaults(use_threads=False)
    p.add_argument(
        '--clean', 
        required=False, 
        dest='clean',
        action='store_true',
        help='Clean the intermediate files.')
    p.set_defaults(clean=False)
    return p.parse_args()


def task1():
    print 'task1'
    with open('task1.txt', 'w') as ofile:
        ofile.write('task1')

def task2(param1):
    print 'task2, param1', param1
    action = 'wc -l task1.txt > task2.txt'
    os.system(action)

def task3():
    print 'task3'
    action = 'bzip2 -k -f task1.txt'
    os.system(action)

def task4():
    print 'task4'
    action = 'rm task1.txt'
    os.system(action)


def task_dependence_example():
    # a task dependence example
    DoitLoader.add_task(
                        targets=['task2.txt'], 
                        file_dep=['task1.txt'],
                        actions=[[task2, {'param1':123}], 'echo HELLO'],
                        name='task2')
 
    DoitLoader.add_task(
                        targets=['task1.txt.bz2'], 
                        file_dep=['task1.txt', 'task2.txt'],
                        actions=task3,
                        name='task3')

    DoitLoader.add_task(
                        targets=[], 
                        file_dep=['task1.txt.bz2'],
                        actions=[task4],
                        name='task4',
                        uptodate=[False])

    DoitLoader.add_task(
                        targets=['task5.txt'], 
                        file_dep=[],
                        actions=[['echo TASK5', {'stdout':'task5.txt'}], 
                                 'echo TASK5_2'],
                        name='task5',
                        uptodate=[False])
 
    DoitLoader.add_task(
                        targets=['task1.txt'], 
                        file_dep=[],
                        actions=[task1],
                        name='task1')


def pipe_task_example():
    # a task built by a pipe of actions
    actions = []
    actions.append([
                    'cat', 
                    {
                     'stdin':'file1.txt'
                    }
                   ])
    actions.append(['grep 12', {'stdout':'file2.txt'}])
    #actions.append(['ls -l'])
    #actions.append('echo abc')
    DoitLoader.add_task(
                          targets=['file2.txt'],
                          file_dep=['file1.txt'],
                          actions=actions,
                          pipe=True)

def pipe_task_example2():
    # a task built by a pipe of actions
    actions = []
    actions.append('ls')
    actions.append([
                    'cat', 
                    {
                     'stdin':'file2.txt',
                     'stdout':subprocess.PIPE
                    }
                   ])
    actions.append(['grep 12', {'stdout':'file3.txt'}])
    #actions.append(['ls -l'])
    actions.append('echo abc')
    DoitLoader.add_task(
                          targets=['file3.txt'],
                          file_dep=['file2.txt'],
                          actions=actions)


def main(args):
    #task_dependence_example()
    #pipe_task_example()
    pipe_task_example2()


    if args.clean:
        doit_args = ['clean']
    else:
        doit_args = ['-n', str(args.nprocs)]
        if args.use_threads:
            doit_args += ['-P', 'thread']

    code = DoitLoader.run(doit_args)
    sys.exit(code)


if __name__ == "__main__":
    args = read_params()
    main(args)
