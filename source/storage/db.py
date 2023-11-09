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
from sqlalchemy.dialects.postgresql import JSONB
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
    #willl need to set up a postgresql db for json
    #charge_model_charges = Column(MutableDict.as_mutable(JSONB), nullable=True)

#class DBChargeModel(DBBase):
#    id = Column(Integer, primary_key=True)
#    conformer_id = Column(Integer, ForeignKey("conformers.id"))
#    charges = Column(MutableDict.as_mutable(JSONB), nullable=True)

class DBMoleculeRecordProp(DBMoleculeRecord):

    conformers = relationship("DBConformerRecordProp")

