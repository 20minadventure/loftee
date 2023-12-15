import subprocess
import pkg_resources


file_path = pkg_resources.resource_filename('preprocessing', 'install_vep.sh')

def main():
    subprocess.run(['bash', file_path], check=True)
