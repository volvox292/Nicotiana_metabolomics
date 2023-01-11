#!/usr/bin/env python3
# coding: utf-8

import subprocess
import argparse
import os
from pathlib import Path
from contextlib import contextmanager
from molmass import Formula
import time

@contextmanager
def timeit(desc):
  start = time.time()
  yield
  print(f'{desc} finished in {int(time.time() - start)} seconds')

@contextmanager
def set_directory(path: Path):
  origin = Path().absolute()
  try:
    os.chdir(path)
    yield
  finally:
    os.chdir(origin)

# def smiles_to_coord(smiles):
#   ifl = 'coord.tmol'
#   subprocess.check_call(['obabel', f'-:{smiles}', '-O', ifl, '--gen3d', 'best'])
#   subprocess.check_call(['obabel', ifl,  '-O', 'coord.pdb'])
#   subprocess.check_call(['mv', ifl, 'coord'])

def inchi_to_coord(inchi):
  ifl = 'coord.tmol'
  subprocess.check_call(['obabel', '-i', 'inchi', f'-:{inchi}', '-O', ifl, '--gen3d', 'best'])
  subprocess.check_call(['obabel', ifl,  '-O', 'coord.pdb'])
  subprocess.check_call(['mv', ifl, 'coord'])

def extract_protomer(fname):
  protomers=[]
  with open(fname, 'r') as f:
    for line in f:
      line=line.strip('\n')
      protomers.append(line)

  x = int(protomers[0]) + 2
  coord_path="coord.xyz"
  with open(coord_path, 'w') as f:
    for i in protomers[0:x]:
      f.write(i)
      f.write('\n')

def xyz_to_coord(xyz):
  subprocess.check_call(['obabel', xyz,  '-O', 'coord.tmol'])
  subprocess.check_call(['mv', 'coord.tmol', 'coord'])

def write_qcxms_config(config):
  with open('qcxms.in', 'w') as f:
    for opt, val in config.items():
      if not val:
        f.write(f'{opt}\n')
      else:
        f.write(f'{opt} {val}\n')

def status_print(s):
  print("\r", s, end="")

def inchi_to_molecular_formula(inchi):
  molecular_formula = inchi.split('/')[1]
  return(molecular_formula)

def elab_predict_from_mass(mass):
  # Prediction made from a set of 4 mass and 4 elab values
  return(round(34.6070 + (0.2359 * mass), 0))

def deduc_elab_from_inchi(inchi):
  # Use molmass to convert molecular formula to molecular mass
  f = Formula(inchi_to_molecular_formula(inchi))
  return(elab_predict_from_mass(f.isotope.mass))

class Commands:
  @staticmethod
  def cid(inchi, outdir, ntraj=10, **kwargs):
    os.makedirs(outdir, exist_ok=True)
    with set_directory(Path(outdir)):
      inchi_to_coord(inchi)
      elab = deduc_elab_from_inchi(inchi)
      write_qcxms_config(
        {
          'xtb2': '',
          'cid': '',
          'elab': str(elab),
          'iatom':'n2',
          'lchamb': '0.25',
          'fullauto': '',
          'ntraj': str(ntraj),
        }
      )

      with timeit(f'QCxMS CID for {inchi}'):
        try:
          status_print('Calculating protomers')
          subprocess.check_output(['crest', '-protonate'])
          extract_protomer("protonated.xyz")
          xyz_to_coord("coord.xyz")
          status_print('Calculating trajectories')
          subprocess.check_output(['qcxms'])
          status_print('Preparing for prod run')
          subprocess.check_output(['qcxms'])
          status_print(f'Running for {ntraj} trajectories')
          subprocess.check_output(['pqcxms'])
          status_print(f'Getting results')
          subprocess.check_output(['getres'])
          status_print(f'PlotMS for extracting spectra')
          subprocess.check_output(['plotms', '-i'])
          print('Finished')
        except subprocess.CalledProcessError as e:
          print(e)
          print(e.output.decode('utf-8'))

def main(args):
  # Run QCxMS
  Commands.cid(args.inchi, args.outdir, ntraj=args.num_traj)

if __name__ == '__main__':

    # Create parser with inchi as input
    parser = argparse.ArgumentParser(description='Running QC-XMS')
    parser.add_argument('-i', '--inchi', help='InChI code (string). Required.', required=True)
    parser.add_argument('-o', '--outdir', help='Output directory (string).', required=True)
    parser.add_argument('-nt', '--num_traj', help='Number of trajectories (int). Default is 5.', default=5)
    args = parser.parse_args()

    #os.system(f'export OMP_NUM_THREADS={threads}')
    main(args)

