# Copyright 2018-2019 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

"""Object Relational Mapping for ITS1 sequence database.

Using SQLalchemy, the Python classes defined here give us a
database schema and the code to import/export the data as
Python objects.
"""

from sqlalchemy import Column, Integer, String, DateTime, ForeignKey
from sqlalchemy import create_engine, Index
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
    md5 = Column(String(32), unique=True)
    date = Column(DateTime)
    notes = Column(String(1000))

    def __repr__(self):
        """Represent a DataSource database entry as a string."""
        return "DataSource(name=%r, ...)" % self.name


class Taxonomy(Base):
    """Database entry for a species' taxonomy entry."""

    __tablename__ = "taxonomy"
    __table_args__ = (
        # Might have old entry with clade A, and new curated entry with clade B
        # (but same genus, species and NCBI taxid)
        Index(
            "taxid_genus_species_clade",
            "ncbi_taxid",
            "genus",
            "species",
            "clade",
            unique=True,
        ),
        Index("taxid_genus_species", "ncbi_taxid", "genus", "species", unique=False),
        Index("genus_species", "genus", "species", unique=False),
    )

    id = Column(Integer, primary_key=True)
    # Using empty string rather than Null (None) for clade, genus, species
    clade = Column(String(10), nullable=False)  # TODO - Integer linked to table?
    ncbi_taxid = Column(Integer)
    genus = Column(String(100), nullable=False)
    species = Column(String(100), nullable=False)  # source may have variant/strain?

    def __repr__(self):
        """Represent a taxonomy database entry as a string."""
        return "Taxonomy(clade=%r, ncbi_taxid=%r, genus=%r, species=%r)" % (
            self.clade,
            self.ncbi_taxid,
            self.genus,
            self.species,
        )


class ITS1(Base):
    """Database entry for a single ITS1 sequence."""

    __tablename__ = "its1_sequence"

    id = Column(Integer, primary_key=True)
    md5 = Column(String(32), unique=True)
    sequence = Column(String(250), unique=True)

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

    its1_id = Column(Integer, ForeignKey("its1_sequence.id"))
    its1 = relationship(ITS1, foreign_keys=[its1_id])

    sequence = Column(String(1000))  # Full sequence, can be longer than ITS1

    # Whatever was recorded in the original data source
    original_taxonomy_id = Column(Integer, ForeignKey("taxonomy.id"))
    original_taxonomy = relationship(Taxonomy, foreign_keys=[original_taxonomy_id])

    # Initially based on the values above, but expect to reclassify some
    current_taxonomy_id = Column(Integer, ForeignKey("taxonomy.id"))
    current_taxonomy = relationship(Taxonomy, foreign_keys=[current_taxonomy_id])

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


def connect_to_db(*args, **kwargs):
    """Create engine and return session make bound to it.

    >>> Session = connect_to_db('sqlite:///:memory:', echo=True)
    >>> session = Session()
    """
    engine = create_engine(*args, **kwargs)
    Base.metadata.create_all(engine)
    return sessionmaker(bind=engine)


if __name__ == "__main__":
    print("Debugging example")
    Session = connect_to_db("sqlite:///:memory:", echo=True)
