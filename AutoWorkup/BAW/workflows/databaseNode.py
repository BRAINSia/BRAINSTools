"""
databaseNode.py
=========================
Description:

Author:

Usage:

"""
import os.path
import sqlite3

from nipype.interfaces.base import (
    TraitedSpec,
    traits,
    File,
    BaseInterface,
    DynamicTraitedSpec,
    BaseInterfaceInputSpec,
)
from nipype.interfaces.io import IOBase


class SQLiteGrabberInputSpec(DynamicTraitedSpec, BaseInterfaceInputSpec):
    """
    This class represents a...

    :param DynamicTraitedSpec:
    :param BaseInterfaceInputSpec:
    """

    database_file = File(exists=True, mandatory=True)
    table_name = traits.Str(mandatory=True)
    columns = traits.List(traits.Str, mandatory=True)
    constraints = traits.List(
        traits.Tuple(traits.Str, traits.Either(traits.Str, traits.List(traits.Str))),
        minlen=1,
        desc="The list of column/value pairs for WHERE creation",
    )
    distinct = traits.Bool(default=False, usedefault=True)
    orderby = traits.List(
        traits.Tuple(traits.Str, traits.Enum(("ASC", "DESC"))),
        desc="List of tuples with column and order direction",
    )


class SQLiteGrabberOutputSpec(TraitedSpec):
    """
    This class represents a...

    :param TraitedSpec:
    """

    results = traits.List(traits.Tuple(), desc="Results of query")
    query = traits.Str(desc="Query sent to database")


class SQLiteGrabber(IOBase):
    """ Very simple frontend for getting values from SQLite database.

        .. warning::

            This is not a thread-safe node.

        .. warning:: Vulnerable to SQL injection attacks - http://en.wikipedia.org/wiki/SQL_injection

        Examples
        --------

        >>> sql = SQLiteGrabber()
        >>> sql.inputs.database_file = 'my_database.db'
        >>> sql.inputs.table_name = 'experiment_results'
        >>> sql.inputs.columns = ['columnA', 'columnB']
        >>> sql.inputs.constraints = [('column1', 'value1'), ('column2', ['value2a', 'value2b'])]
        >>> sql.query
        SELECT columnA, columnB FROM experiment_results WHERE column1='value1' AND column2 in ('value2a', 'value2b')
        >>> sql.run() # doctest: +SKIP

    """

    input_spec = SQLiteGrabberInputSpec
    output_spec = SQLiteGrabberOutputSpec
    _always_run = True

    def __init__(self, **inputs):
        """This function...

        :param **inpusts:
        """
        self._query = ""
        super(SQLiteGrabber, self).__init__(**inputs)

    @property
    def query(self):
        """ `query` plus any columns and constraints ()
        validates arguments and generates command line"""
        self._check_mandatory_inputs()
        return self._gen_query()

    def _check_mandatory_inputs(self):
        """
        This function will...
        """
        assert os.path.splitext(self.inputs.database_file)[-1] == ".db"
        super(SQLiteGrabber, self)._check_mandatory_inputs()

    def _list_outputs(self):
        """Execute this module.
        """
        outputs = self.output_spec().get()
        outputs["query"] = self.query
        outputs["results"] = self._execute_query()
        return outputs

    def _execute_query(self):
        """
        This funciton...
        """
        conn = sqlite3.connect(self.inputs.database_file, check_same_thread=False)
        c = conn.cursor()
        query = self.query
        try:
            c.execute(query)
            retval = c.fetchall()
        except sqlite3.Error as err:
            # err.message = err.message + " -> " + query
            raise err
        finally:  # Clean up...
            try:
                c.close()
                conn.close()
            except:
                raise  # TODO: Give warning of premature disconnect
        return retval

    def _gen_query(self):
        """
        This function...

        :return:
        """
        # TODO: write SQL query to prevent injection attacks
        if self.inputs.distinct:
            _select = "SELECT DISTINCT"
        else:
            _select = "SELECT"
        query = " ".join(
            [
                _select,
                ", ".join(["{0}".format(c) for c in self.inputs.columns]),
                "FROM {table}".format(table=self.inputs.table_name),
            ]
        )
        if self.inputs.constraints:
            query += " WHERE"
            for key, value in self.inputs.constraints:
                if isinstance(value, str) or isinstance(value, str):
                    query += " {column}='{value}'".format(column=key, value=value)
                elif isinstance(value, list):
                    query += " {column} IN (".format(column=key)
                    query += ", ".join(["'{value}'".format(value=v) for v in value])
                    query += ")"
                query += " AND"
            query = query[:-4]  # Remove last unnecessary " AND"
        if self.inputs.orderby:
            # Remove last unnecessary ","
            query = " ".join(
                [query, "ORDER BY"]
                + [
                    " ".join([column, order.upper() + ","])
                    for column, order in self.inputs.orderby
                ]
            )[:-1]
        return query


def open_subject_database(
    ExperimentBaseDirectoryCache, single_subject, mountPrefix, subject_data_file
):
    """
    This function does...

    :param ExperimentBaseDirectoryCache:
    :param single_subject:
    :param mountPrefix:
    :param subject_data_file:
    :return:
    """
    import os.path
    import SessionDB

    subjectDatabaseFile = os.path.join(
        ExperimentBaseDirectoryCache, "InternalWorkflowSubjectDB.db"
    )
    ## TODO:  Only make DB if db is older than subject_data_file.
    if (not os.path.exists(subjectDatabaseFile)) or (
        os.path.getmtime(subjectDatabaseFile) < os.path.getmtime(subject_data_file)
    ):
        ExperimentDatabase = SessionDB.SessionDB(subjectDatabaseFile, single_subject)
        ExperimentDatabase.make_new_db(subject_data_file, mountPrefix)
        ExperimentDatabase = None
        ExperimentDatabase = SessionDB.SessionDB(subjectDatabaseFile, single_subject)
    else:
        print(
            (
                "Single_subject {0}: Using cached database, {1}".format(
                    single_subject, subjectDatabaseFile
                )
            )
        )
        ExperimentDatabase = SessionDB.SessionDB(subjectDatabaseFile, single_subject)
    # print "ENTIRE DB for {_subjid}: ".format(_subjid=ExperimentDatabase.get_subject_filter())
    # print "^^^^^^^^^^^^^"
    # for row in ExperimentDatabase.get_everything():
    #    print row
    # print "^^^^^^^^^^^^^"
    return ExperimentDatabase


def get_all_scans(cache, subject, prefix, dbfile, session):
    """
    This function...

    :param cache:
    :param subject:
    :param prefix:
    :param dbfile:
    :param session:
    :return:
    """
    pass


def make_database_node(
    cache, dbfile, table_name="MasterDB", columns=["*"], constraints=[]
):
    """
    This function...

    :param cache:
    :param dbfile:
    :param table_name:
    :param columns:
    :param constraints:
    :return:
    """
    import os.path
    import nipype.pipeline.engine as pe  # pypeline engine
    from .databaseNode import SQLiteGrabber  # open_subject_database

    node = pe.Node(interface=SQLiteGrabber(), name="99_OpenSubjectDatabase")
    node.inputs.database_file = os.path.join(cache, dbfile)
    node.inputs.table_name = table_name
    node.inputs.columns = columns
    node.inputs.constraints = constraints
    return node


def session_constraint(session):
    """
    This function...

    :param session:
    :return:
    """
    return [("session", session)]


def files_constraint(session, types):
    """
    This function...

    :param session:
    :param types:
    :return:
    """
    return [("session", session), ("type", types)]
