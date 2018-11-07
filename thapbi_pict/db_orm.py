"""Object Relational Mapping for ITS1 sequence database.

Using SQLalchemy, the Python classes defined here give us a
database schema and the code to import/export the data as
Python objects.
"""

from sqlalchemy import Column, Integer, String, DateTime, ForeignKey
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
    its1_md5 = Column(String(32), ForeignKey("its1_sequence.md5"))
    sequence = Column(String(1000))  # Full sequence, can be longer than ITS1

    # Whatever was recorded in the original data source
    original_taxid = Column(Integer)  # NCBI taxid
    original_genus = Column(String(100))
    original_species = Column(String(100))

    # Initially the same values as above, but expect to reclassify some entries
    current_taxid = Column(Integer)  # NCBI taxid
    current_genus = Column(String(100))
    current_species = Column(String(100))

    # TODO - Supplement this with an NCBI taxonomy table?
    # i.e. looking up species name and genus via the taxid
    # For now, storing genus etc locally will help with simple filtering

    date_added = Column(DateTime)
    date_modified = Column(DateTime)

    # Implement as an enum?
    # Unknown, Public DB, Direct Sequencing, Genome Sequence, Metabarcoding
    seq_strategy = Column(Integer)

    # Implement as an enum?
    # Unknown, Sanger, Illumina, Roche 454, PacBio, Oxford Nanopore, etc
    seq_platform = Column(Integer)

    # Start with something like -1 = bad, +1 = good?
    curated_trust = Column(Integer)


Base.metadata.create_all(engine)
