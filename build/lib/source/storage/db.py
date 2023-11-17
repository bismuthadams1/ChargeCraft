from openff.recharge.esp.storage.db import DBMoleculeRecord, DBGridSettings, DBPCMSettings, DBESPSettings, DBConformerRecord, DBBase

from sqlalchemy import (
    Boolean,
    Column,
    ForeignKey,
    Integer,
    PickleType,
    String,
    UniqueConstraint
)
from sqlalchemy.orm import Query, Session, relationship
#implemenet postgresql in the future 
#from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.ext.mutable import MutableDict
from sqlalchemy.ext.declarative import declarative_base





class DBConformerRecordProp(DBConformerRecord):
    mulliken_charges = Column(PickleType, nullable=False)
    lowdin_charges = Column(PickleType, nullable=False)
    mbis_charges  = Column(PickleType, nullable=False)
    dipole = Column(PickleType, nullable=False)
    quadropole =  Column(PickleType, nullable=False)
    mbis_dipole = Column(PickleType, nullable=False)
    mbis_quadropole = Column(PickleType, nullable=False)
    mbis_octopole = Column(PickleType, nullable=False)
    energy = Column(PickleType, nullable=False)
    charge_model_charges = Column(String, nullable=True) # Change JSONB to String

class DBMoleculeRecordProp(DBMoleculeRecord):

    conformers = relationship("DBConformerRecordProp")

