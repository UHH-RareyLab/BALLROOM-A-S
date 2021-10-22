import pandas as pd
from sqlalchemy import create_engine, text

from src.utils import definitions


def query_chembl(query):
    # Erstellen Sie die SQLAlchemy-Engine
    engine = create_engine(
        f"postgresql+psycopg2://{definitions.chembl_username}:{definitions.chembl_password}@whakaari/{definitions.chembl_databasename}"
    )

    with engine.connect() as connection:
        result_df = pd.read_sql(text(query), con=connection)
    return result_df


def query_activitydb(query):
    engine = create_engine(
        f"postgresql+psycopg2://{definitions.activitydb_username}:{definitions.activitydb_password}@whakaari/{definitions.activitydb_databasename}"
    )

    with engine.connect() as connection:
        result_df = pd.read_sql(text(query), con=connection)
    return result_df


def query_or_readin_data(df_path, sql_query, recalc, chembl=True):
    if df_path.exists() and not recalc:
        return pd.read_csv(df_path, index_col=0)
    else:
        if chembl:
            current_data = query_chembl(sql_query)
        else:
            current_data = query_activitydb(sql_query)
        current_data.to_csv(df_path)
        return current_data


def query_activities(chembl_mol, chembl_target):
    """
    This function queries the activities for a specific chembl mol and chembl target in our version of chembl
    :param chembl_mol: Chembl mol identifier, starts with CHEMBL and continues with a number
    :param chembl_target: Chembl target identifier, starts with CHEMBL and continues with a number
    :return:
    """
    sql = f"""SELECT 
    a.activity_id AS activity_id,
    a.standard_value,
    a.standard_units,
    t.chembl_id AS target_id,
    t.target_type,
    m.chembl_id AS molecule_id,
    m.pref_name,
    assays.chembl_id
    FROM 
        activities a
    JOIN molecule_dictionary m ON a.molregno= m.molregno
    JOIN assays ON a.assay_id= assays.assay_id
    JOIN target_dictionary t ON assays.tid = t.tid
    WHERE 
        m.chembl_id = '{chembl_mol}' AND 
        t.chembl_id = '{chembl_target}';"""
    dataframe = query_chembl(sql)
    return dataframe
