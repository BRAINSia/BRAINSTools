from __future__ import print_function


def run_workflow(workflow, plugin='Linear', plugin_args={}):
    """
    Run workflow object and catch traceback for printing to stdout
    """
    import traceback
    import sys

    print("Running workflow...")
    try:
        workflow.run(plugin=plugin, plugin_args=plugin_args)
    except:
        print("=+-+" * 25)
        print("Error: Exception while running subjects")
        traceback.print_exc(file=sys.stdout)
        return False
    return True


def print_workflow(workflow, plugin, dotfilename='workflow', graph2use='hierarchical'):
    """
    HINT: graph2use values: ['orig', 'flat', 'hierarchical', 'exec']
    'hierarchical' is the only one that DOES NOT require pygraphviz
    """
    assert plugin in ['Linear', 'MultiProc'], "'plugin' must be in ['Linear', 'MultiProc'] to print workflow"
    dotfilename = '_'.join([dotfilename, graph2use])
    print("Writing graph to filename {0}".format(dotfilename))
    try:
        workflow.write_graph(dotfilename=dotfilename, graph2use=graph2use)
        return True
    except:
        raise
    return False
