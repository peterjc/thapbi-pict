"""Object Relational Mapping for ITS1 sequence database.

Using SQLalchemy, the Python classes defined here give us a
database schema and the code to import/export the data as
Python objects.
"""

from sqlalchemy import Column, Integer, String, DateTime, ForeignKey
from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, sessionmaker


Base = declarative_base()


class DataSource(Base):
    """Database entry for a data source (NCBI, Legacy, etc).

    Each accession is expected to be unique within a data source.
    """

    __tablename__ = "data_source"

    id = Column(Integer, primary_key=True)
    name = Column(String(100))  # e.g. NCBI, Legacy v0.005
    uri = Column(String(200))  # e.g. traceable filename or URL
    date = Column(DateTime)
    notes = Column(String(1000))

    def __repr__(self):
        """Represent a DataSource database entry as a string."""
        return "DataSource(name=%r, ...)" % self.name


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

    id = Column(Integer, primary_key=True)

    source_accession = Column(String)  # Hopefully unique within source_id
    source_id = Column(Integer, ForeignKey("data_source.id"))
    source = relationship(DataSource)

    its1_md5 = Column(String(32), ForeignKey("its1_sequence.md5"))
    sequence = Column(String(1000))  # Full sequence, can be longer than ITS1

    # Whatever was recorded in the original data source
    original_clade = Column(String(10))  # TODO - Integer linked to table?
    original_taxid = Column(Integer)  # NCBI taxid
    original_genus = Column(String(100))
    original_species = Column(String(100))  # source may have variant/strain?

    # Initially based on the values above, but expect to reclassify some
    current_clade = Column(String(10))
    current_taxid = Column(Integer)  # NCBI taxid, perhaps as species only?
    current_genus = Column(String(100))
    current_species = Column(String(100))

    # TODO - Supplement this with an NCBI taxonomy table?
    # i.e. looking up species name and genus via the taxid
    # For now, storing genus etc locally will help with simple filtering

    # date_added = Column(DateTime) -> see data_source.date
    date_modified = Column(DateTime)

    # Implement as an enum?
    # Unknown, Public DB, Direct Sequencing, Genome Sequence, Metabarcoding
    seq_strategy = Column(Integer)

    # Implement as an enum?
    # Unknown, Sanger, Illumina, Roche 454, PacBio, Oxford Nanopore, etc
    seq_platform = Column(Integer)

    # Start with something like -1 = bad, +1 = good?
    curated_trust = Column(Integer)

    its1_seq = relationship("ITS1", back_populates="entries")


ITS1.entries = relationship("SequenceSource",
                            order_by=SequenceSource.source_accession,
                            back_populates="its1_seq")


def connect_to_db(*args, **kwargs):
    """Create engine and return session make bound to it.

    >>> Session = connect_to_db('sqlite:///:memory:', echo=True)
    >>> session = Session()
    """
    engine = create_engine(*args, **kwargs)
    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    return Session


if __name__ == "__main__":
    print("Debugging example")
    Session = connect_to_db('sqlite:///:memory:', echo=True)
