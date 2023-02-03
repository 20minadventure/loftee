import pytest
import collections
import dxpy

from analysis.utils.dxpathlib import PathDx


@pytest.fixture
def test_database_name():
    return 'test_database'


@pytest.fixture
def test_database_id():
    return 'database-GPJ7YJjJ18XYJy56jZgF7FZ6'


@pytest.fixture
def local_non_empty_dir_path():
    return PathDx() / 'tests'


@pytest.fixture
def dnax_non_empty_dir_path(test_database_id):
    return PathDx(database_id=test_database_id)


@pytest.fixture
def dnax_empty_table_path(test_database_id):
    return PathDx('empty_table.hl', database_id=test_database_id)


def test_docker_paths():
    assert PathDx().rstr == 'file:///opt/notebooks/ifpan-gosborcz-drugs'
    assert PathDx('example_table.ht').rstr == 'file:///opt/notebooks/ifpan-gosborcz-drugs/example_table.ht'
    assert PathDx('/example_table.ht').rstr == 'file:///example_table.ht'


def test_dnax_paths(test_database_name, test_database_id):
    example_table_real_path = f'dnax://{test_database_id}/example_table.ht'
    
    p1 = PathDx('example_table.ht', database=test_database_name)
    p2 = PathDx('/example_table.ht', database=test_database_name)
    p3 = PathDx(database=test_database_name) / 'example_table.ht'
    p4 = PathDx('example_table.ht', database_id=test_database_id)

    assert p1.rstr == example_table_real_path
    assert p2.rstr == example_table_real_path
    assert p3.rstr == example_table_real_path
    assert p4.rstr == example_table_real_path


def test_non_existing_database_name():
    with pytest.raises(ValueError):
        PathDx('example_table.ht', database='i_dont_exist')
    with pytest.raises(ValueError):
        str(PathDx('example_table.ht', database='i_dont_exist'))


def test_database_specification(test_database_name, test_database_id):
    with pytest.raises(ValueError):
        PathDx(database=test_database_name, database_id=test_database_id)

    path_by_name = PathDx(database=test_database_name).rstr
    path_by_id = PathDx(database=test_database_id).rstr
    assert path_by_name == path_by_id


def test_iterdir(local_non_empty_dir_path, dnax_non_empty_dir_path):
    assert isinstance(
        local_non_empty_dir_path.iterdir(),
        collections.Iterator
    )
    assert isinstance(
        dnax_non_empty_dir_path.iterdir(),
        collections.Iterator
    )
    
    assert isinstance(
        next(local_non_empty_dir_path.iterdir()),
        PathDx
    )
    assert isinstance(
        next(dnax_non_empty_dir_path.iterdir()),
        PathDx
    )


def test_listdir(
        local_non_empty_dir_path, dnax_non_empty_dir_path, test_database_id):
    assert len(local_non_empty_dir_path.listdir()) > 0
    assert len(dnax_non_empty_dir_path.listdir()) > 0
      
    table_p = PathDx('empty_table.ht', database_id=test_database_id)
    assert len(table_p.listdir()) == 7

    one_dir_p = PathDx('empty_table.ht/index', database_id=test_database_id)
    assert len(one_dir_p.listdir()) == 1
    
    file_p = PathDx('empty_table.ht/README.txt', database_id=test_database_id)
    with pytest.raises(NotADirectoryError):
        file_p.listdir()
    
    empty_dir_p = PathDx(
        'empty_table.ht/references',
        database_id=test_database_id
    )
    with pytest.raises(Warning):
        empty_dir_p.listdir()
    
    non_exising_p = PathDx(
        'empty_table.ht/i_dont_exist',
        database_id=test_database_id
    )
    with pytest.raises(Warning):
        non_exising_p.listdir()
