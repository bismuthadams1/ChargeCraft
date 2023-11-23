from openff.recharge.esp.storage.db import (DBMoleculeRecord, 
                                            DBGridSettings, 
                                            DBPCMSettings, 
                                            DBESPSettings, 
                                            DBConformerRecord, 
                                            DBBase,
                                            _UniqueMixin,
                                            _DB_FLOAT_PRECISION,
                                            _float_to_db_int,
                                            _db_int_to_float)

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

from source.storage.ddx_storage import DDXSettings


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

    ddx_settings = relationship("DBDDXSettings", uselist = False)
    ddx_settings_id = Column(Integer, ForeignKey("ddx_settings.id"), nullable = True)


class DBMoleculeRecordProp(DBMoleculeRecord):

    conformers = relationship("DBConformerRecordProp")

class DBDDXSettings(_UniqueMixin, DBBase):
    __tablename__ = "ddx_settings"

    id = Column(Integer, primary_key=True, index=True)
 
    ddx_model = Column(String(6), nullable = False)
    solvent = Column(String(20), nullable = True)
    epsilon = Column(Integer, nullable = True)
    radii_set = Column(String(5), nullable = False)

    @classmethod
    def _hash(cls, instance: DDXSettings) -> int:
        return hash(
            (
                instance.ddx_model,
                instance.solvent,
                _float_to_db_int(instance.epsilon),
                instance.radii_set
            )
        )

    @classmethod
    def _query(cls, db: Session, instance: DDXSettings) -> Query:
        epsilon = _float_to_db_int(instance.solvent)

        return (
            db.query(DBDDXSettings)
            .filter(DBDDXSettings.ddx_model == instance.ddx_model)
            .filter(DBDDXSettings.solvent == instance.solvent)
            .filter(DBDDXSettings.epsilon == instance.epsilon)
            .filter(DBDDXSettings.radii_set == instance.radii_set)
        )

    @classmethod
    def _instance_to_db(cls, instance: DDXSettings) -> "DBDDXSettings":
        return DBDDXSettings(
            ddx_model=instance.ddx_model,
            solvent=instance.solvent,
            epsilon=instance.epsilon,
            radii_set=instance.radii_set
        )

    @classmethod
    def db_to_instance(cls, db_instance: "DBDDXSettings") -> DDXSettings:
        # noinspection PyTypeChecker
        return PCMSettings(
            ddx_model=db_instance.ddx_model,
            solvent=db_instance.solvent,
            epsilon=_float_to_db_int(db_instance.epsilon),
            radii_set=db_instance.radii_set
        )   









