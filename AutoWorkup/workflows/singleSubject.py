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
    start_time, subject, master_config = args
    assert 'baseline' in master_config['components'] or 'longitudinal' in master_config[
        'components'], "Baseline or Longitudinal is not in WORKFLOW_COMPONENTS!"
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
    from utils import run_workflow, print_workflow

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
    # To avoid a "sqlite3.ProgrammingError: Base Cursor.__init__ not called" error
    #    using multiprocessing.map_async(), instantiate database here
    database = OpenSubjectDatabase(
        master_config['cachedir'], [subject], master_config['prefix'], master_config['dbfile'])
    # print database.getAllSessions()
    database.open_connection()

    sessions = database.getSessionsFromSubject(subject)
    # print "These are the sessions: ", sessions
    if 'baseline' in master_config['components']:
        current_phase = 'baseline'
    elif 'longitudinal' in master_config['components']:
        current_phase = 'longitudinal'
    from longitudinal import create_longitudinal as create_wkfl

    # TODO (future): Replace with iterable inputSpec node and add Function node for getAllFiles()
    for session in sessions:
        project = database.getProjFromSession(session)
        # Long node names make graphs a pain to read/print
        pname = "{0}_{1}".format(session, current_phase)
        # pname = GenerateWFName(project, subject, session, current_phase)
        print "Building session pipeline for {0}".format(session)
        inputsSpec[session] = pe.Node(name='inputspec_{0}'.format(session),
                                      interface=IdentityInterface(fields=['T1s', 'T2s', 'PDs', 'FLs', 'OTs']))
        inputsSpec[session].inputs.T1s = database.getFilenamesByScantype(
            session, ['T1-15', 'T1-30'])
        inputsSpec[session].inputs.T2s = database.getFilenamesByScantype(
            session, ['T2-15', 'T2-30'])
        inputsSpec[session].inputs.PDs = database.getFilenamesByScantype(
            session, ['PD-15', 'PD-30'])
        inputsSpec[session].inputs.FLs = database.getFilenamesByScantype(
            session, ['FL-15', 'FL-30'])
        inputsSpec[session].inputs.OTs = database.getFilenamesByScantype(
            session, ['OTHER-15', 'OTHER-30'])

        sessionWorkflow[session] = create_wkfl(project, subject, session, master_config,
                                               interpMode='Linear', pipeline_name=pname)

        subjectWorkflow.connect(
            [(inputsSpec[session], sessionWorkflow[session], [('T1s', 'inputspec.T1s'),
                                                              ('T2s',
                                                               'inputspec.T2s'),
                                                              ('PDs',
                                                               'inputspec.PDs'),
                                                              ('FLs',
                                                               'inputspec.FLs'),
                                                              ('OTs',
                                                               'inputspec.OTHERs'),
                                                              ]),
             (atlasNode, sessionWorkflow[
                 session], [('template_landmarks_50Lmks_fcsv',
                             'inputspec.atlasLandmarkFilename'),
                            ('template_weights_50Lmks_wts',
                             'inputspec.atlasWeightFilename'),
                            ('LLSModel_50Lmks_hdf5',
                             'inputspec.LLSModel'),
                            ('T1_50Lmks_mdl', 'inputspec.inputTemplateModel')]),
             ])
        if 'segmentation' in master_config['components']:
            from WorkupT1T2BRAINSCut import GenerateWFName
            try:
                bCutInputName = ".".join(
                    ['segmentation', GenerateWFName(project, subject, session, 'Segmentation'), 'inputspec'])
            except:
                print project, subject, session
                raise
            subjectWorkflow.connect([(atlasNode, sessionWorkflow[session],
                                      [('hncma-atlas', 'segmentation.inputspec.hncma-atlas'),
                                       ('template_t1', 'segmentation.inputspec.template_t1'),
                                       ('template_t1', bCutInputName + '.template_t1'),
                                       ('rho', bCutInputName + '.rho'),
                                       ('phi', bCutInputName + '.phi'),
                                       ('theta', bCutInputName + '.theta'),
                                       ('l_caudate_ProbabilityMap',
                                        bCutInputName + '.l_caudate_ProbabilityMap'),
                                       ('r_caudate_ProbabilityMap',
                                        bCutInputName + '.r_caudate_ProbabilityMap'),
                                       ('l_hippocampus_ProbabilityMap',
                                        bCutInputName + '.l_hippocampus_ProbabilityMap'),
                                       ('r_hippocampus_ProbabilityMap',
                                        bCutInputName + '.r_hippocampus_ProbabilityMap'),
                                       ('l_putamen_ProbabilityMap',
                                        bCutInputName + '.l_putamen_ProbabilityMap'),
                                       ('r_putamen_ProbabilityMap',
                                        bCutInputName + '.r_putamen_ProbabilityMap'),
                                       ('l_thalamus_ProbabilityMap',
                                        bCutInputName + '.l_thalamus_ProbabilityMap'),
                                       ('r_thalamus_ProbabilityMap',
                                        bCutInputName + '.r_thalamus_ProbabilityMap'),
                                       ('l_accumben_ProbabilityMap',
                                        bCutInputName + '.l_accumben_ProbabilityMap'),
                                       ('r_accumben_ProbabilityMap',
                                        bCutInputName + '.r_accumben_ProbabilityMap'),
                                       ('l_globus_ProbabilityMap',
                                        bCutInputName + '.l_globus_ProbabilityMap'),
                                       ('r_globus_ProbabilityMap',
                                        bCutInputName + '.r_globus_ProbabilityMap'),
                                       ('trainModelFile_txtD0060NT0060_gz',
                                        bCutInputName + '.trainModelFile_txtD0060NT0060_gz')])])
            subjectWorkflow.connect(
                [(atlasNode, sessionWorkflow[session], [('template_t1', 'inputspec.template_t1'),
                                                        ('ExtendedAtlasDefinition_xml',
                                                         'inputspec.atlasDefinition')]),
                 ])
