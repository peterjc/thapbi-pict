"""Object Relational Mapping for ITS1 sequence database.

Using SQLalchemy, the Python classes defined here give us a
database schema and the code to import/export the data as
Python objects.
"""

from sqlalchemy import Column, Integer, String
from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base


engine = create_engine('sqlite:///:memory:', echo=True)
Base = declarative_base()


class ITS1(Base):
    """Database entry for a single ITS1 sequence."""

    __tablename__ = "its1_sequence"

    md5 = Column(String(32), primary_key=True)
    sequence = Column(String(250))

    def __repr__(self):
        """Represent an ITS1 database entry as a string."""
        return "ITS1(md5=%r, sequence=%r)" % (self.md5, self.sequence)


class SequenceSource(Base):
    """Database entry for source of an ITS1 sequence entry."""

    __tablename__ = "its1_source"

    accession = Column(String, primary_key=True)
    its1_md5 = Column(String(32))
    sequence = Column(String(1000))  # Full sequence, can be longer than ITS1
    taxid = Column(Integer)  # NCBI taxid
    species = Column(String(100))


Base.metadata.create_all(engine)
