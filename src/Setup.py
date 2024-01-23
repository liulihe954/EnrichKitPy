import subprocess
import os
from pathlib import Path

class Setup:
    '''
    A class for setting up the base environment, including downloading SQLite and more.
    '''
    def __init__(self, base, zenodo_get, zenodo_record_id, cloud_store='zenodo'):
        self.base = Path(base)
        self.zenodo_get = zenodo_get
        self.zenodo_record_id = zenodo_record_id
        self.cloud_store = cloud_store

    def download(self):
        dbcache_path = self.base / '.dbcache'
        db_path = dbcache_path / 'EnrichKitDB.sqlite'

        if not dbcache_path.exists() or not db_path.exists():
            if not dbcache_path.exists():
                dbcache_path.mkdir(parents=True, exist_ok=True)
            if self.cloud_store == 'zenodo':
                subprocess.run([self.zenodo_get, self.zenodo_record_id], check=True)
                print('Download completed!')
        else:
            print('Core DB/SQLite exists!')

