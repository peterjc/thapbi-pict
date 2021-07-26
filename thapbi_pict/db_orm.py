# Copyright 2018-2021 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.
"""Object Relational Mapping for marker sequence database.

Using SQLalchemy, the Python classes defined here give us a
database schema and the code to import/export the data as
Python objects.
"""
from sqlalchemy import Column
from sqlalchemy import create_engine
from sqlalchemy import ForeignKey
from sqlalchemy import Index
from sqlalchemy import Integer
from sqlalchemy import String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy.orm import sessionmaker


Base = declarative_base()


class DataSource(Base):
    """Database entry for a data source (NCBI, curated, etc).

    Each accession is expected to be unique within a data source.
    """

    __tablename__ = "data_source"

    id = Column(Integer, primary_key=True)
    name = Column(String(100))  # e.g. NCBI, Legacy v0.005
    uri = Column(String(200))  # e.g. traceable filename or URL
    md5 = Column(String(32), unique=True)
    notes = Column(String(1000))

    def __repr__(self):
        """Represent a DataSource database entry as a string."""
        return f"DataSource(name={self.name!r}, ...)"


class Taxonomy(Base):
    """Database entry for a species' taxonomy entry."""

    __tablename__ = "taxonomy"
    __table_args__ = (
        Index("taxid_genus_species", "ncbi_taxid", "genus", "species", unique=True),
        Index("genus_species", "genus", "species", unique=True),
    )

    id = Column(Integer, primary_key=True)
    # Using empty string rather than Null (None) for genus, species
    ncbi_taxid = Column(Integer)
    genus = Column(String(100), nullable=False)
    species = Column(String(100), nullable=False)  # source may have variant/strain?

    def __repr__(self):
        """Represent a taxonomy database entry as a string."""
        return (
            f"Taxonomy(ncbi_taxid={self.ncbi_taxid!r},"
            f" genus={self.genus!r}, species={self.species!r})"
        )


class Synonym(Base):
    """Database entry for a synonym of a taxonomy entry.

    In addition to direct synonyms, includes the names and synonyms of any
    child nodes of the species (e.g. variants, strains, etc).
    """

    __tablename__ = "synonym"

    id = Column(Integer, primary_key=True)
    taxonomy_id = Column(Integer, ForeignKey("taxonomy.id"))
    name = Column(String(100), nullable=False, unique=True)  # genus and species in one

    def __repr__(self):
        """Represent a synonym database entry as a string."""
        return f"Synonym(name={self.name})"


class MarkerDef(Base):
    """Database entry for a marker listing primers and amplicon length limits."""

    __tablename__ = "marker_definition"

    id = Column(Integer, primary_key=True)
    name = Column(String(32), unique=True)
    left_primer = Column(String(100))
    right_primer = Column(String(100))
    min_length = Column(Integer)
    max_length = Column(Integer)


class MarkerSeq(Base):
    """Database entry for a single marker reference sequence."""

    __tablename__ = "marker_sequence"

    id = Column(Integer, primary_key=True)
    md5 = Column(String(32), unique=True)
    sequence = Column(String(250), unique=True)

    def __repr__(self):
        """Represent a marker database reference sequence as a string."""
        return f"MarkerSeq(md5={self.md5!r}, sequence={self.sequence!r})"


class SeqSource(Base):
    """Database entry for source of a marker sequence entry."""

    __tablename__ = "sequence_source"

    id = Column(Integer, primary_key=True)

    source_accession = Column(String)  # Hopefully unique within source_id
    source_id = Column(Integer, ForeignKey("data_source.id"))
    source = relationship(DataSource)

    marker_seq_id = Column(Integer, ForeignKey("marker_sequence.id"))
    marker_seq = relationship(MarkerSeq, foreign_keys=[marker_seq_id])

    marker_definition_id = Column(Integer, ForeignKey("marker_definition.id"))
    marker_definition = relationship(MarkerDef, foreign_keys=[marker_definition_id])

    taxonomy_id = Column(Integer, ForeignKey("taxonomy.id"))
    taxonomy = relationship(Taxonomy, foreign_keys=[taxonomy_id])


def connect_to_db(*args, **kwargs):
    """Create engine and return session bound to it.

    >>> Session = connect_to_db('sqlite:///:memory:', echo=True)
    >>> session = Session()
    """
    engine = create_engine(*args, **kwargs)
    Base.metadata.create_all(engine)
    return sessionmaker(bind=engine)


if __name__ == "__main__":
    print("Debugging example")
    Session = connect_to_db("sqlite:///:memory:", echo=True)
