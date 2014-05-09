def run_workflow(workflow, plugin='Linear', plugin_args={}):
    """
    Run workflow object and catch traceback for printing to stdout
    """
    import traceback
    import sys

    print "Running workflow..."
    try:
        workflow.run(plugin_args=plugin_args)
        # if not plugin in ['Linear', 'MultiProc']:
        #     workflow.run(plugin_args=plugin_args)
        # else:
        #     workflow.run()
    except:
        print "=+-+" * 25
        print "Error: Exception while running subjects"
        traceback.print_exc(file=sys.stdout)
        return False
    return True

def print_workflow(workflow, plugin, dotfilename='workflow', graph2use='hierarchial'):
    if plugin in ['Linear', 'MultiProc']:
        print "Writing graph to filename {0}".format(dotfilename)
        try:
            workflow.write_graph(dotfilename=dotfilename, graph2use=graph2use)
            return True
        except:
            raise
    return False
