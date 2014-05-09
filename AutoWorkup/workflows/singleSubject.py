def RunSubjectWorkflow(args):
    """
                           .-----------.
                       --- | Session 1 | ---> /project/subjectA/session1/phase/
                     /     *-----------*
    .-----------.   /
    | Subject A | <
    *-----------*   \
                     \     .-----------.
                       --- | Session 2 | ---> /project/subjectA/session2/phase/
                           *-----------*
    **** Replaces WorkflowT1T2.py ****
    """
    database, start_time, subject, master_config = args
    assert 'baseline' in master_config['components'] or 'longitudinal' in master_config['components'], "Baseline or Longitudinal is not in WORKFLOW_COMPONENTS!"
    # HACK:
    #    To avoid a "sqlite3.ProgrammingError: Base Cursor.__init__ not called" error
    #    using multiprocessing.map_async(), re-instantiate database
    # database.__init__(defaultDBName=database.dbName, subject_list=database.subjectList)
    #
    # END HACK
    import time

    from nipype import config, logging
    config.update_config(master_config)  # Set universal pipeline options
    assert config.get('execution', 'plugin') == master_config['execution']['plugin']
    # DEBUG
    # config.enable_debug_mode()
    # config.set('execution', 'stop_on_first_rerun', 'true')
    # END DEBUG
    logging.update_logging(config)

    import nipype.pipeline.engine as pe
    import nipype.interfaces.base as nbase
    import nipype.interfaces.io as nio
    from nipype.interfaces.utility import IdentityInterface, Function
    import traits

    from baw_exp import OpenSubjectDatabase
    from SessionDB import SessionDB
    from PipeLineFunctionHelpers import convertToList
    from atlasNode import MakeAtlasNode
    from utilities.misc import GenerateSubjectOutputPattern as outputPattern
    from utilities.misc import GenerateWFName

    while time.time() < start_time:
        time.sleep(start_time - time.time() + 1)
        print "Delaying start for {subject}".format(subject=subject)
    print("===================== SUBJECT: {0} ===========================".format(subject))

    subjectWorkflow = pe.Workflow(name="BAW_StandardWorkup_subject_{0}".format(subject))
    subjectWorkflow.base_dir = config.get('logging', 'log_directory')
    # subjectWorkflow.config['execution']['plugin'] = 'Linear'  # Hardcodeded in WorkupT1T2.py - why?
    # DEBUG
    # subjectWorkflow.config['execution']['stop_on_first_rerun'] = 'true'
    # END DEBUG
    atlasNode = MakeAtlasNode(master_config['atlascache'], 'BAtlas')

    sessionWorkflow = dict()
    inputsSpec = dict()
    sessions = database.getSessionsFromSubject(subject)
    # print "These are the sessions: ", sessions
    if 'baseline' in master_config['components']:
        current_phase = 'baseline'
        from baseline import create_baseline as create_wkfl
    elif 'longitudinal' in master_config['components']:
        current_phase = 'longitudinal'
        from longitudinal import create_longitudial as create_wkfl

    for session in sessions:  # TODO (future): Replace with iterable inputSpec node and add Function node for getAllFiles()
        project = database.getProjFromSession(session)
        pname = "{0}_{1}".format(session, current_phase)  # Long node names make graphs a pain to read/print
        # pname = GenerateWFName(project, subject, session, current_phase)
        print "Building session pipeline for {0}".format(session)
        inputsSpec[session] = pe.Node(name='inputspec_{0}'.format(session),
                                      interface=IdentityInterface(fields=['T1s', 'T2s', 'PDs', 'FLs', 'OTs']))
        inputsSpec[session].inputs.T1s = database.getFilenamesByScantype(session, ['T1-15', 'T1-30'])
        inputsSpec[session].inputs.T2s = database.getFilenamesByScantype(session, ['T2-15', 'T2-30'])
        inputsSpec[session].inputs.PDs = database.getFilenamesByScantype(session, ['PD-15', 'PD-30'])
        inputsSpec[session].inputs.FLs = database.getFilenamesByScantype(session, ['FL-15', 'FL-30'])
        inputsSpec[session].inputs.OTs = database.getFilenamesByScantype(session, ['OTHER-15', 'OTHER-30'])

        sessionWorkflow[session] = create_wkfl(project, subject, session, master_config,
                                               interpMode='Linear', pipeline_name=pname)

        subjectWorkflow.connect([(inputsSpec[session], sessionWorkflow[session], [('T1s', 'inputspec.T1s'),
                                                                                  ('T2s', 'inputspec.T2s'),
                                                                                  ('PDs', 'inputspec.PDs'),
                                                                                  ('FLs', 'inputspec.FLs'),
                                                                                  ('OTs', 'inputspec.OTHERs'),
                                                                                  ]),
                                 (atlasNode, sessionWorkflow[session], [('template_landmarks_50Lmks_fcsv',
                                                                         'inputspec.atlasLandmarkFilename'),
                                                                        ('template_weights_50Lmks_wts',
                                                                         'inputspec.atlasWeightFilename'),
                                                                        ('LLSModel_50Lmks_hdf5', 'inputspec.LLSModel'),
                                                                        ('T1_50Lmks_mdl', 'inputspec.inputTemplateModel')]),
                                ])
        if current_phase == 'baseline':
            subjectWorkflow.connect([(atlasNode, sessionWorkflow[session], [('template_t1', 'inputspec.template_t1'),
                                                                            ('ExtendedAtlasDefinition_xml',
                                                                             'inputspec.atlasDefinition')]),
                                 ])
        else:
            assert current_phase == 'longitudinal', "Phase value is unknown: {0}".format(current_phase)

    from utils import run_workflow, print_workflow
    if False:
        print_workflow(template, plugin=master_config['execution']['plugin'], dotfilename='template')
    return run_workflow(template, plugin=master_config['execution']['plugin'], plugin_args=master_config['plugin_args'])
