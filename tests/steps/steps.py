import os
from behave import *
import shlex
import subprocess
from contextlib import contextmanager


root_dir = os.path.normpath(
    os.path.join(os.path.dirname(__file__), '..', '..'))
print('ROOT DIR', root_dir)


@given('build executable {filename} exists')
def step_impl(context, filename):
    """ Check existance in root directory """
    filename = os.path.join(root_dir, filename)
    if os.path.exists(filename):
        context.executable = filename
    else:
        raise FileNotFoundError(filename)


@given('the file {filename} exists')
def step_impl(context, filename):
    """ Check existance in root directory """
    filename = os.path.join(root_dir, filename)
    if os.path.exists(filename):
        pass
    else:
        raise FileNotFoundError(filename)


@when('I run the executable with arguments {arguments}')
def step_impl(context, arguments):
    args = shlex.split(arguments)
    command = [context.executable]
    command.extend(args)
    with pushd(root_dir):
        print('Running command', command, 'in', os.getcwd())
        completed_process = subprocess.run(command)
    context.return_code = completed_process.returncode


@then('the command is successful')
def step_impl(context):
    if context.return_code == 0:
        pass
    else:
        raise AssertionError(
            'Exit code was {}, expected 0'.format(context.return_code))


@then('the file {filename} is created')
def step_impl(context, filename):
    filename = os.path.join(root_dir, filename)
    if os.path.exists(filename):
        pass
    else:
        raise FileNotFoundError(filename)


@contextmanager
def pushd(directory):
    old_dir = os.getcwd()
    os.chdir(directory)
    yield
    os.chdir(old_dir)
