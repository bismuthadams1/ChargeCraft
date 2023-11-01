from openff.recharge.esp.storage.db import DBMoleculeRecord, DBGridSettings, DBPCMSettings, DBESPSettings, DBConformerRecord
from sqlalchemy import (
    Boolean,
    Column,
    ForeignKey,
    Integer,
    PickleType,
    String,
    UniqueConstraint,
)
from sqlalchemy.orm import Query, Session, relationship



class DBConformerRecordProp(DBConformerRecord):
    mulliken_charges = Column(PickleType, nullable=False)
    lowdin_charges = Column(PickleType, nullable=False)
    mbis_charges  = Column(PickleType, nullable=False)
    dipole = Column(PickleType, nullable=False)
    quadropole =  Column(PickleType, nullable=False)


class DBMoleculeRecordProp(DBMoleculeRecord):

    conformers = relationship("DBConformerRecordProp")

