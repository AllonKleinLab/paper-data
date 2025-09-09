import os, subprocess
import itertools
import operator
from collections import defaultdict, OrderedDict
import errno

try:
   import cPickle as pickle
except:
   import pickle

import numpy as np
import re
import shutil
import gzip
from itertools import product, combinations
import time
import yaml
import tempfile
import string
from contextlib import contextmanager

# -----------------------
# Helper functions
# -----------------------

def string_hamming_distance(str1, str2):
    """Fast hamming distance over 2 strings known to be of same length."""
    return sum(itertools.imap(operator.ne, str1, str2))

def rev_comp(seq):
    """Reverse complement of DNA sequence."""
    tbl = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    return ''.join(tbl[s] for s in seq[::-1])

def to_fastq(name, seq, qual):
    """Return string that can be written to fastQ file."""
    return '@'+name+'\n'+seq+'\n+\n'+qual+'\n'

def seq_neighborhood(seq, n_subs=1):
    """Generate all sequences within n_subs substitutions of input sequence."""
    for positions in combinations(range(len(seq)), n_subs):
        for subs in product(*("ATGCN",)*n_subs):
            seq_copy = list(seq)
            for p, s in zip(positions, subs):
                seq_copy[p] = s
            yield ''.join(seq_copy)

def build_barcode_neighborhoods(barcode_file, expect_reverse_complement=True):
    """Build barcode neighborhoods for error correction."""
    clean_mapping = dict()
    mapping1 = defaultdict(set)
    mapping2 = defaultdict(set)
    
    with open(barcode_file, 'rU') as f:
        for line in f:
            barcode = line.rstrip()
            if expect_reverse_complement:
                barcode = rev_comp(line.rstrip())

            clean_mapping[barcode] = barcode

            for n in seq_neighborhood(barcode, 1):
                mapping1[n].add(barcode)
            
            for n in seq_neighborhood(barcode, 2):
                mapping2[n].add(barcode)   
    
    for k, v in mapping1.items():
        if k not in clean_mapping:
            if len(v) == 1:
                clean_mapping[k] = list(v)[0]
    
    for k, v in mapping2.items():
        if k not in clean_mapping:
            if len(v) == 1:
                clean_mapping[k] = list(v)[0]
    
    return clean_mapping

def check_dir(path):
    """Create directory if it doesn't exist."""
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def print_to_stderr(msg, newline=True):
    """Write message to stderr."""
    import sys
    sys.stderr.write(str(msg))
    if newline:
        sys.stderr.write('\n')

def worker_filter(iterable, worker_index, total_workers):
    """Filter items for specific worker."""
    return (p for i,p in enumerate(iterable) if (i-worker_index)%total_workers==0)

class FIFO():
    """Context manager for a named pipe."""
    def __init__(self, filename="", suffix="", prefix="tmp_fifo_dir", dir=None):
        if filename:
            self.filename = filename
        else:
            self.tmpdir = tempfile.mkdtemp(suffix=suffix, prefix=prefix, dir=dir)
            self.filename = os.path.join(self.tmpdir, 'fifo')

    def __enter__(self):
        if os.path.exists(self.filename):
            os.unlink(self.filename)
        os.mkfifo(self.filename)
        return self

    def __exit__(self, type, value, traceback):
        os.remove(self.filename)
        if hasattr(self, 'tmpdir'):
            shutil.rmtree(self.tmpdir)

# -----------------------
# Core objects
# -----------------------

class IndropsProject():
    def __init__(self, project_yaml_file_handle):
        self.yaml = yaml.safe_load(project_yaml_file_handle)
        self.name = self.yaml['project_name']
        self.project_dir = self.yaml['project_dir']
        self.libraries = OrderedDict()
        self.runs = OrderedDict()

        for run in self.yaml['sequencing_runs']:
            version = run['version']
            
            if version != 'v3' and version != 'v3-miseq':
                continue  # Only handle V3 for filtering
                
            filtered_filename = '{library_name}_{run_name}_{library_index}'
            if 'split_affixes' in run:
                filtered_filename += '_{split_affix}'
                split_affixes = run['split_affixes']
            else:
                split_affixes = ['']

            filtered_filename += '.fastq'

            if 'libraries' in run:
                run_libraries = run['libraries']
            elif 'library_name' in run:
                run_libraries = [{'library_name' : run['library_name'], 'library_prefix':''}]
            else:
                raise Exception('No library name or libraries specified.')

            for affix in split_affixes:
                filtered_part_filename = filtered_filename.format(run_name=run['name'], split_affix=affix,
                    library_name='{library_name}', library_index='{library_index}')
                part_filename = os.path.join(self.project_dir, '{library_name}', 'filtered_parts', filtered_part_filename)

                input_filename = os.path.join(run['dir'], run['fastq_path'].format(split_affix=affix, read='{read}'))
                part = V3Demultiplexer(run['libraries'], project=self, part_filename=part_filename, 
                                     input_filename=input_filename, run_name=run['name'], part_name=affix,
                                     run_version_details=run['version'])

                if run['name'] not in self.runs:
                    self.runs[run['name']] = []
                self.runs[run['name']].append(part)

                for lib in run_libraries:
                    lib_name = lib['library_name']
                    lib_index = lib['library_index']
                    if lib_name not in self.libraries:
                        self.libraries[lib_name] = IndropsLibrary(name=lib_name, project=self, version=run['version'])
                    self.libraries[lib_name].parts.append(part.libraries[lib_index])

    @property
    def paths(self):
        if not hasattr(self, '_paths'):
            script_dir = os.path.dirname(os.path.realpath(__file__))
            with open(os.path.join(script_dir, 'default_parameters.yaml'), 'r') as f:
                paths = yaml.safe_load(f)['paths']
            paths.update(self.yaml['paths'])

            paths['python'] = os.path.join(paths['python_dir'], 'python')
            paths['java'] = os.path.join(paths['java_dir'], 'java')
            paths['trimmomatic_jar'] = os.path.join(script_dir, 'bins', 'trimmomatic-0.33.jar')

            self._paths = type('Paths_anonymous_object',(object,),paths)()
            self._paths.trim_polyA_and_filter_low_complexity_reads_py = os.path.join(script_dir, 'trim_polyA_and_filter_low_complexity_reads.py')
            self._paths.count_barcode_distribution_py = os.path.join(script_dir, 'count_barcode_distribution.py')
            self._paths.gel_barcode1_list = os.path.join(script_dir, 'ref/barcode_lists/gel_barcode1_list.txt')
            self._paths.gel_barcode2_list = os.path.join(script_dir, 'ref/barcode_lists/gel_barcode2_list.txt')
        return self._paths

    @property
    def parameters(self):
        if not hasattr(self, '_parameters'):
            with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'default_parameters.yaml'), 'r') as f:
                self._parameters = yaml.safe_load(f)['parameters']
            if 'parameters' in self.yaml:
                for k, d in self.yaml['parameters'].items():
                    self._parameters[k].update(d)
        return self._parameters

    @property
    def gel_barcode1_revcomp_list_neighborhood(self):
        if not hasattr(self, '_gel_barcode1_revcomp_list_neighborhood'):
            self._gel_barcode1_revcomp_list_neighborhood = build_barcode_neighborhoods(self.paths.gel_barcode1_list, False)
        return self._gel_barcode1_revcomp_list_neighborhood
  
    @property
    def gel_barcode2_list_neighborhood(self):
        if not hasattr(self, '_gel_barcode2_list_neighborhood'):
            self._gel_barcode2_list_neighborhood = build_barcode_neighborhoods(self.paths.gel_barcode2_list, False)
        return self._gel_barcode2_list_neighborhood


class IndropsLibrary():
    def __init__(self, name='', project=None, version=''):
        self.project = project
        self.name = name
        self.parts = []
        self.version = version

        self.paths = {}
        for lib_dir in ['filtered_parts']:
            dir_path = os.path.join(self.project.project_dir, self.name, lib_dir)
            check_dir(dir_path)
            self.paths[lib_dir] = dir_path
        self.paths = type('Paths_anonymous_object',(object,),self.paths)()


class LibrarySequencingPart():
    def __init__(self, filtered_fastq_filename=None, project=None, run_name='', library_name='', part_name=''):
        self.project = project
        self.run_name = run_name
        self.part_name = part_name
        self.library_name = library_name
        self.filtered_fastq_filename = filtered_fastq_filename
        self.barcode_counts_pickle_filename = filtered_fastq_filename + '.counts.pickle'
        self.filtering_metrics_filename = '.'.join(filtered_fastq_filename.split('.')[:-1]) + 'metrics.yaml'
        self.filtering_statistics_counter = defaultdict(int)

    def contains_library_in_query(self, query_libraries):
        return self.library_name in query_libraries

    @contextmanager
    def trimmomatic_and_low_complexity_filter_process(self):
        """Start trimmomatic and complexity filter processes."""
        filtered_dir = os.path.dirname(self.filtered_fastq_filename)
        
        with FIFO(dir=filtered_dir) as fifo2, open(self.filtered_fastq_filename, 'w') as filtered_fastq_file, open(self.filtered_fastq_filename+'.counts.pickle', 'w') as filtered_index_file:
            
            low_complexity_filter_cmd = [self.project.paths.python, self.project.paths.trim_polyA_and_filter_low_complexity_reads_py,
                '-input', fifo2.filename, 
                '--min-post-trim-length', self.project.parameters['trimmomatic_arguments']['MINLEN'],
                '--max-low-complexity-fraction', str(self.project.parameters['low_complexity_filter_arguments']['max_low_complexity_fraction']),
                ]
            counter_cmd = [self.project.paths.python, self.project.paths.count_barcode_distribution_py]

            p2 = subprocess.Popen(low_complexity_filter_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            p3 = subprocess.Popen(counter_cmd, stdin=p2.stdout, stdout=filtered_fastq_file, stderr=filtered_index_file)

            with FIFO(dir=filtered_dir) as fifo1:
                trimmomatic_cmd = [self.project.paths.java, '-Xmx500m', '-jar', self.project.paths.trimmomatic_jar,
                        'SE', '-threads', "1", '-phred33', fifo1.filename, fifo2.filename]
                for arg in self.project.parameters['trimmomatic_arguments']['argument_order']:
                    val = self.project.parameters['trimmomatic_arguments'][arg]
                    trimmomatic_cmd.append('%s:%s' % (arg, val))

                p1 = subprocess.Popen(trimmomatic_cmd, stderr=subprocess.PIPE)

                fifo1_filehandle = open(fifo1.filename, 'w')
                yield fifo1_filehandle
                fifo1_filehandle.close()
                
                trimmomatic_stderr = p1.stderr.read().splitlines()
                if trimmomatic_stderr[2] != 'TrimmomaticSE: Completed successfully':
                    raise Exception('Trimmomatic did not complete succesfully')
                trimmomatic_metrics = trimmomatic_stderr[1].split() 
                trimmomatic_metrics = {'input' : trimmomatic_metrics[2], 'output': trimmomatic_metrics[4], 'dropped': trimmomatic_metrics[7]}
                p1.wait()

            complexity_filter_metrics = pickle.load(p2.stderr)
            p2.wait()
            p3.wait()

        filtering_metrics = {
            'read_structure' : dict(self.filtering_statistics_counter),
            'trimmomatic' : trimmomatic_metrics,
            'complexity_filter': complexity_filter_metrics,
        }
        with open(self.filtering_metrics_filename, 'w') as f:
            yaml.dump(dict(filtering_metrics), f, default_flow_style=False)


class V3Demultiplexer():
    def __init__(self, library_indices, project=None, part_filename="", input_filename="", run_name="", part_name="", run_version_details="v3"):
        self.run_version_details = run_version_details
        self.input_filename = input_filename
        self.project = project
        self.run_name = run_name
        self.part_name = part_name
        self.libraries = {}
        for lib in library_indices:
            lib_index = lib['library_index']
            lib_name = lib['library_name']
            library_part_filename = part_filename.format(library_name=lib_name, library_index=lib_index)
            self.libraries[lib_index] = LibrarySequencingPart(filtered_fastq_filename=library_part_filename, 
                                                              project=project, run_name=run_name, 
                                                              library_name=lib_name, part_name=part_name)

    def _weave_fastqs(self, fastqs):
        """Read and weave multiple FASTQ files together."""
        last_extension = [fn.split('.')[-1] for fn in fastqs]
        if all(ext == 'gz' for ext in last_extension):
            processes = [subprocess.Popen("gzip --stdout -d %s" % (fn), shell=True, stdout=subprocess.PIPE) for fn in fastqs]
            streams = [r.stdout for r in processes]
        elif all(ext == 'bz2' for ext in last_extension):
            processes = [subprocess.Popen("bzcat %s" % (fn), shell=True, stdout=subprocess.PIPE) for fn in fastqs]
            streams = [r.stdout for r in processes]
        elif all(ext == 'fastq' for ext in last_extension):
            streams = [open(fn, 'r') for fn in fastqs]
        else:
            raise("ERROR: Different files are compressed differently. Check input.")

        while True:
            names = [next(s)[:-1].split()[0] for s in streams]
            seqs = [next(s)[:-1] for s in streams]
            blanks = [next(s)[:-1]  for s in streams]
            quals = [next(s)[:-1]  for s in streams]
            assert all(name==names[0] for name in names)
            yield names[0], seqs, quals

        for s in streams:
            s.close()

    def _process_reads(self, name, seqs, quals, valid_bc1s={}, valid_bc2s={}, valid_libs={}):
        """Process reads for V3 format."""
        r1, r2, r3, r4 = seqs

        if r3 in valid_libs:
            lib_index = valid_libs[r3]
        else:
            return False, r3, 'Invalid_library_index'
        
        orig_bcd1 = r2[:6] + r4[21:27]
        if orig_bcd1 in valid_bc1s:
            bc1 = valid_bc1s[orig_bcd1]
        else:
            return False, lib_index, 'Invalid_BC1'

        orig_bc2 = r4[:6] + r4[10:16]
        umi = r4[27:37]

        if orig_bc2 in valid_bc2s:
            bc2 = valid_bc2s[orig_bc2]
        else:
            return False, lib_index, 'Invalid_BC2'

        if 'N' in umi:
            return False, lib_index, 'UMI_contains_N'

        final_bc = '%s-%s' % (bc1, bc2)
        return True, lib_index, (final_bc, umi)

    def filter_and_count_reads(self):
        """Main filtering function for V3 format."""
        # Prepare error corrected index sets
        self.sequence_to_index_mapping = {}
        libs = self.libraries.keys()
        self.sequence_to_index_mapping = dict(zip(libs, libs))
        index_neighborhoods = [set(seq_neighborhood(lib, 1)) for lib in libs]
        for lib, clibs in zip(libs, index_neighborhoods):
            for clib in clibs:
                if sum(clib in hood for hood in index_neighborhoods)==1:
                    self.sequence_to_index_mapping[clib] = lib

        # Prepare error corrected barcode sets
        error_corrected_barcodes = self.project.gel_barcode2_list_neighborhood
        error_corrected_rev_compl_barcodes = self.project.gel_barcode1_revcomp_list_neighborhood

        # Open up context managers
        manager_order = []
        trim_processes = {}
        trim_processes_managers = {}
        orig_fastq_handles = {}
        barcode_fastq_handles = {}  

        for lib in self.libraries.keys():
            manager_order.append(lib)
            trim_processes_managers[lib] = self.libraries[lib].trimmomatic_and_low_complexity_filter_process()
            trim_processes[lib] = trim_processes_managers[lib].__enter__()

            orig_fastq = self.libraries[lib].filtered_fastq_filename
            if not orig_fastq.endswith('.gz'):
                orig_fastq += '.gz'
            orig_fastq_handles[lib] = gzip.open(orig_fastq, 'wt')
            
            # Open barcode+UMI FastQ for writing (gzipped)
            barcode_fastq = self.libraries[lib].filtered_fastq_filename.replace('.fastq', '_barcodes.fastq')
            if not barcode_fastq.endswith('.gz'):
                barcode_fastq += '.gz'
            barcode_fastq_handles[lib] = gzip.open(barcode_fastq, 'wt')

        overall_filtering_statistics = defaultdict(int)
        input_fastqs = []
        for r in ['R1', 'R2', 'R3', 'R4']:
            input_fastqs.append(self.input_filename.format(read=r))

        # Setup logging
        last_ping = time.time()
        ping_every_n_reads = 1000000
        ping_header = "{0:>12}{1:>16}{2:>12}{3:>10}{4:>10}{5:>10}{6:>10}   |" + ''.join("{%d:>12.10}"%i for i in range(7,7+len(manager_order)))
        ping_header = ping_header.format("Total Reads", "", "Valid Reads", "No index", "No BC1", "No BC2", "No UMI", *[self.libraries[k].library_name for k in manager_order])
        ping_template = "{total:12d}    {rate:5.1f} sec/M {Valid:12.1%}{Invalid_library_index:10.1%}{Invalid_BC1:10.1%}{Invalid_BC2:10.1%}{UMI_contains_N:10.1%}   |{"+":>12.1%}{".join(manager_order)+":>12.1%}"
        
        def print_ping_to_log(last_ping):
            sec_per_mil = (time.time() - last_ping)/(float(ping_every_n_reads)/10**6) if last_ping else 0
            total = overall_filtering_statistics['Total']
            if total > 0:
                ping_format_data = {k: float(overall_filtering_statistics[k])/total for k in ['Valid', 'Invalid_library_index', 'Invalid_BC1',  'Invalid_BC2', 'UMI_contains_N']}
                if overall_filtering_statistics['Valid'] > 0:
                    ping_format_data.update({k: float(self.libraries[k].filtering_statistics_counter['Valid'])/overall_filtering_statistics['Valid'] for k in manager_order})
                else:
                    ping_format_data.update({k: 0.0 for k in manager_order})
                print_to_stderr(ping_template.format(total=total, rate=sec_per_mil, **ping_format_data))

        len_bc1, len_bc2, len_umi = 12, 12, 10  # Barcode and UMI lengths

        print_to_stderr('Filtering %s, file %s' % (self.run_name, self.input_filename))

        for r_name, seqs, quals in self._weave_fastqs(input_fastqs):
            # Python 3 compatibility
            seqs = [s.decode('utf-8') if isinstance(s, bytes) else s for s in seqs]

            keep, lib_index, result = self._process_reads(r_name, seqs, quals,
                                                        error_corrected_rev_compl_barcodes, 
                                                        error_corrected_barcodes, 
                                                        self.sequence_to_index_mapping)
            if keep:
                bc, umi = result
                bio_read = seqs[0]
                bio_qual = quals[0]
                read_id_extra = r_name[1:]  # get the rest of the Illumina read header line after '@'

                # Compose headers for both outputs
                orig_header = "{}:{} {}".format(bc, umi, read_id_extra)
                barcode_header = "{}:{} {}".format(bc, umi, read_id_extra)

                orig_fastq_handles[lib_index].write(to_fastq(orig_header, bio_read, bio_qual))

                # Barcode+UMI output
                try:
                    bc1, bc2 = bc.split('-')
                except Exception:
                    print("Barcode splitting failed for:", bc)
                    bc1, bc2 = bc, ''
                barcode_seq = bc1 + bc2 + umi
                barcode_qual = bio_qual[:len_bc1+len_bc2+len_umi]

                barcode_fastq_handles[lib_index].write(to_fastq(barcode_header, barcode_seq, barcode_qual))

                self.libraries[lib_index].filtering_statistics_counter['Valid'] += 1
                self.libraries[lib_index].filtering_statistics_counter['Total'] += 1
                overall_filtering_statistics['Valid'] += 1

            else:
                if result != 'Invalid_library_index':
                    self.libraries[lib_index].filtering_statistics_counter[result] += 1
                    self.libraries[lib_index].filtering_statistics_counter['Total'] += 1
                overall_filtering_statistics[result] += 1

            overall_filtering_statistics['Total'] += 1

            # Print logging header every 10M reads
            if overall_filtering_statistics['Total']%(ping_every_n_reads*10)==1:
                print_to_stderr(ping_header)
            
            # Print progress every 1M reads
            if overall_filtering_statistics['Total']%ping_every_n_reads == 0:
                print_ping_to_log(last_ping)
                last_ping = time.time()

        # Final log output
        print_ping_to_log(False)

        # Close up the context managers and barcode FASTQ handles
        for lib in manager_order[::-1]:
            trim_processes_managers[lib].__exit__(None, None, None)
            orig_fastq_handles[lib].close()
            barcode_fastq_handles[lib].close()

    def contains_library_in_query(self, query_libraries):
        for lib in self.libraries.values():
            if lib.contains_library_in_query(query_libraries):
                return True
        return False


if __name__=="__main__":
    import sys, argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('project', type=argparse.FileType('r'), help='Project YAML File.')
    parser.add_argument('-l', '--libraries', type=str, help='Library name(s) to work on.', nargs='?', default='')
    parser.add_argument('-r', '--runs', type=str, help='Run name(s) to work on.', nargs='?', default='')
    parser.add_argument('command', type=str, choices=['filter'])
    parser.add_argument('--total-workers', type=int, help='Total workers working together.', default=1)
    parser.add_argument('--worker-index', type=int, help='Index of current worker (starting from 0).', default=0)

    args = parser.parse_args()
    project = IndropsProject(args.project)

    target_libraries = []
    if args.libraries:
        for lib in args.libraries.split(','):
            assert lib in project.libraries
            if lib not in target_libraries:
                target_libraries.append(lib)
    else:
        target_libraries = project.libraries.keys()
    lib_query = set(target_libraries)

    target_runs = []
    if args.runs:
        for run in args.runs.split(','):
            assert run in project.runs
            target_runs.append(run)
    else:
        target_runs = project.runs.keys()

    if args.command == 'filter':
        target_run_parts = []
        for run in target_runs:
            target_run_parts += [part for part in project.runs[run] if part.contains_library_in_query(lib_query)]

        for part in worker_filter(target_run_parts, args.worker_index, args.total_workers):
            print_to_stderr('Filtering run "%s", part "%s"' % (part.run_name, part.part_name))
            part.filter_and_count_reads()
