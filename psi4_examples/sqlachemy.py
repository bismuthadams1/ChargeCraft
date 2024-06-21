from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from chargecraft.storage.storage import MoleculePropRecord, MoleculePropStore
from chargecraft.storage.db import (
    DB_VERSION,
    DBBase,
    DBGridSettings,
    DBGeneralProvenance,
    DBInformation,
    DBMoleculePropRecord,
    DBPCMSettings,
    DBSoftwareProvenance,
    DBConformerPropRecord,
    DBESPSettings,
    DBDDXSettings
)


database_url = 'sqlite:////Users/localadmin/Documents/projects/ChargeCraft/examples/ddxpcmtest.db'
engine = create_engine(database_url, echo=False)

Session = sessionmaker(bind = engine)
session = Session()
session.query(DBMoleculePropRecord)