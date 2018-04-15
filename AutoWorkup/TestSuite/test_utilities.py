from builtins import range

from past.builtins import execfile

from AutoWorkup import utilities


def configure_env_test():
    config_env = os.path.join(os.path.dirname(utilities.__file__), 'configure_env.py')
    for p in range(10):
        file_template = '/my/test/path/{0}'.format(p)

    exec(compile(open(config_env).read(), config_env, 'exec'), dict(__file__=__file__,
                              append_os_path=['/my/test/path/1:/my/test/path/2'],
                              append_os_path=['/my/test/path/3:/my/test/path/4']))
    assert os.environ['PATH'].split(os.pathsep)[0]
