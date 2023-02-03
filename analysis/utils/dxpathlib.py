from pathlib import PosixPath, PurePath

import dxpy


class PathDx(PosixPath):
    DRV_DNAX = 'dnax://'
    DRV_LOCAL = 'file://'

    def __new__(cls, *args, database=None, database_id=None):
        if database is None and database_id is None:
            c = super().__new__(cls, *args)
        elif database is not None and database_id is not None:
            raise ValueError('Both database and database_id are specified.')
        else:
            if database is not None:
                database_id = PathDx.find_database(database)['id']
            c = super().__new__(cls, '/', *args)
            c._drv = f"{PathDx.DRV_DNAX}{database_id}"     
        return c

    @classmethod
    def find_database(cls, db_ref=None):
        databases_dx = dxpy.api.system_find_databases()
        assert databases_dx['next'] is None
        # TODO: handle 'next' request
        databases_dx = databases_dx['results']
        for db in databases_dx:
            database_id = db['id']
            database_desc = dxpy.api.database_describe(database_id)
            database_name = database_desc['name']
            if database_id == db_ref or database_name == db_ref:
                db.update(database_desc)
                return db
        else:
            msg = f"Database doesn't exist. " \
            f"Create with: SC.sql(f\"CREATE DATABASE " \
            f"IF NOT EXISTS {db_ref} LOCATION 'dnax://'\")"
            raise ValueError(msg)

    @property
    def rstr(self):
        database = None
        if self._drv:
            p = str(self)
        else:
            p =  f'{PathDx.DRV_LOCAL}{str(self.resolve())}'
        return p

    def iterdir(self):
        if self._drv:
            database_id = self._drv[len(PathDx.DRV_DNAX):]
            folder = str(PathDx(*self.parts[1:])) if len(self.parts) > 1 else '/'
            files_response = dxpy.api.database_list_folder(
                database_id,
                input_params = {
                    'folder': folder,
                    'includeHidden': True,
                }
            )
            files = [file_d['path'] for file_d in files_response['results']]
            if not files_response['results']:
                raise Warning('Directory is empty OR path does not exist.')

            for result in files_response['results']:
                file_path_split = result['path'].split(database_id)
                file_path = file_path_split[-1]
                if len(file_path_split) == 1:
                    raise NotADirectoryError(f'Not a directory: {file_path}')
                elif len(file_path_split) == 2:
                    yield PathDx(file_path, database_id=database_id)
                else:
                    raise NotImplementedError
        else:
            yield from super().iterdir()

    def listdir(self):
        return list(self.iterdir())
