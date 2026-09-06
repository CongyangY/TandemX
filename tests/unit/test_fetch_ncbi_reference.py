import pytest

from benchmarks.scripts.fetch_ncbi_reference import assembly_parent, resolve_directory


def test_resolves_only_exact_version_from_source_listing():
    accession = 'GCA_034140825.1'
    parent = assembly_parent(accession)
    assert parent.endswith('/GCA/034/140/825/')
    listing = '<a href="../">Parent</a><a href="GCA_034140825.2_later/">new</a><a href="GCA_034140825.1_AGIS1.0/">exact</a>'
    assert resolve_directory(parent, accession, listing) == parent+'GCA_034140825.1_AGIS1.0/'
    with pytest.raises(ValueError, match='unique'):
        resolve_directory(parent, accession, listing+'<a href="GCA_034140825.1_ambiguous/">x</a>')
    with pytest.raises(ValueError, match='unique'):
        resolve_directory(parent, accession, '<a href="https://example.org/GCA_034140825.1_evil/">x</a>')


@pytest.mark.parametrize('value', ['GCA_034140825', 'GCF_034140825.1', 'GCA_034140825.0', '../GCA_034140825.1'])
def test_rejects_unversioned_or_other_assembly_identifiers(value):
    with pytest.raises(ValueError, match='versioned'):
        assembly_parent(value)
