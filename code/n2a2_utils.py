"""
n2a2_utils.py

Contains functions for processing data from the Non-natural Aptamer Array (N2A2)
[REFERENCE]. Generate sequence-intensity data based on provided input data


@author Leighton Wan, Stanford University, 2019
Based on MATLAB code from Peter Mage

"""

import sys
import os
import glob
import numpy as np
import time
import gzip
from pathlib import Path

### GENERAL FILE FUNCTIONS ###

def rename_cycle_directories(run_path,cycle_start,cycle_names):
    """
    Renames the cif directories
    
    Input:
    ---------------------
    run_path: file path to run data
    cycle_start: number of cycle to start
    cycle_names: names to be appended on to the existing directories
    
    Outputs:
    ---------------------
    Renames and prints out status
    
    """
    
    # Prepare format
    cifs_dir=os.path.join(run_path,'cifs')
    cif_dir_format='C{}.1'

    # Find all old directories
    rename_old_new={}
    all_dir_items=os.listdir(cifs_dir)
    # Loop through all the directores
    n_cycles=len(cycle_names)
    for i,cycle_n in enumerate(range(cycle_start,cycle_start+n_cycles)):
        cif_dir_name=cif_dir_format.format(cycle_n)
        old_dir_found=False
        for item in all_dir_items:
            if cif_dir_name in item:
                old_dir=os.path.join(cifs_dir,item)
                new_dir=os.path.join(cifs_dir,cif_dir_name+'_'+cycle_names[i])
                rename_old_new[old_dir]=new_dir
                all_dir_items.remove(item)
                old_dir_found=True
                break
        if not old_dir_found:
            print('No directory found for cycle {}'.format(cycle_n))

    # Rename directories
    for old_dir_i in rename_old_new:
        new_dir_i=rename_old_new[old_dir_i]
        try:
            os.rename(old_dir_i,new_dir_i)
            print('Renamed to: {}'.format(new_dir_i))
        except:
            print('Rename failed for cycle {}'.format(cycle_n))

### LOCS ###

def gen_locs_filename(run_path,tile_num):
    """
    Generates the filename for the loc file
    
    Input:
    ---------------------
    run_path: file path to run data
    tile_num: the tile to use
    
    Outputs:
    ---------------------
    locs_file: path and file name
    
    """
    loc_file='s_1_{}.locs'.format(tile_num)
    locs_file=os.path.join(run_path,'locs',loc_file)
    return locs_file


def parse_locs_file(locs_file):
    """
    Parse locs file into x and y positions
    
    Input:
    ---------------------
    locs_file: String with file and path for locs file
    
    Outputs:
    ---------------------
    xpos,ypos: Positions from locs file as float32
    
    """
    with open(locs_file, 'rb') as file:
        # Remove first twelve bytes
        locs_extra=np.fromfile(file,dtype=np.float32,count=3)
        # Parse the rest
        locs_points=np.fromfile(file,dtype=np.float32)

    # Extract as points
    return locs_points[::2],locs_points[1::2]

### FASTQ ###

def fastq_separate_extract(run_path,fastq_name,verbose=True):
    """
    Parses and separates fastq files into tiles, keeping only sequence and position data
    
    Input:
    ---------------------
    run_path: file path to run data
    fastq_name: name of fastq to process
    verbose: whether to print status updates
    
    Outputs:
    ---------------------
    Writes out files
    
    """
    fastq_path=os.path.join(run_path,'fastq')
    fastq_name=fastq_name+'_L001_R1_001'
    fastq_file=os.path.join(fastq_path,fastq_name+'.fastq.gz')

    # Make directory to store files
    tile_data_dir=os.path.join(fastq_path,fastq_name+'_tile_data')
    if not os.path.exists(tile_data_dir):
        os.makedirs(tile_data_dir)

    with gzip.GzipFile(fastq_file, 'rb') as gzfile:
        tile_current,file_write='',None
        
        if verbose:
            print('Working on FASTQ file: {}'.format(fastq_name))
            t_tile_start=time.time()
        
        for line in gzfile:
            # Extract and decode data
            seq_id=line.rstrip().decode('utf-8')
            seq=gzfile.readline().rstrip().decode('utf-8')
            q_score_id=gzfile.readline().rstrip().decode('utf-8')
            q_score=gzfile.readline().rstrip().decode('utf-8')

            # Extract sequence info (tile,x,y)
            seq_info=seq_id.split(' ')[0]
            seq_parsed=seq_info.split(':')
            tile,seq_x,seq_y=seq_parsed[4:7]

            # Export to appropriate file
            if tile==tile_current:
                file_write.write('{},{},{}\n'.format(seq,seq_x,seq_y))
            else:
                # Close previous
                if file_write is not None:
                    file_write.close()

                # Update tile and open new file
                tile_current=tile
                tile_file='{}_{}.csv'.format(fastq_name,tile_current)
                tile_file_full=os.path.join(tile_data_dir,tile_file)
                file_write=open(tile_file_full,'a+')
                
                # Print update
                if verbose:
                    sys.stdout.write('\r - Current Tile: {}'.format(tile_current))
                    sys.stdout.flush()
                    
        if verbose:
            sys.stdout.write('\n - Completed after {:.2f} seconds\n'.format(time.time()-t_tile_start))
            sys.stdout.flush()

def gen_fastq_filename(run_path,fastq_name,tile_num):
    """
    Generates the filename for the fastq file
    
    Input:
    ---------------------
    run_path: file path to run data
    fastq_name: name of fastq
    tile_num: the tile to use
    
    Outputs:
    ---------------------
    tile_file: path and file name for tile (str)
    
    """
    fastq_dir=fastq_name+'_L001_R1_001_tile_data'
    tile_filename='{}_L001_R1_001_{}.csv'.format(fastq_name,tile_num)
    tile_file=os.path.join(run_path,'fastq',fastq_dir,tile_filename)
    return tile_file


def parse_tile_fastq(tile_file):
    """
    Parse tile fastq into sequence and x,y positions
    
    Input:
    ---------------------
    tile_file: String with file and path for tile fastq file
    
    Outputs:
    ---------------------
    seqs: List of strings for sequences
    xpos,ypos: Array of positions from fastq file as int8
    
    """
    seqs,xpos,ypos=[],[],[]
    with open(tile_file,'r') as file:
        err_count=0
        for line in file:
            seq,x,y=line.rstrip().split(',')
            seqs.append(seq)
            try:
                xpos.append(int(x))
                ypos.append(int(y))
            except:
                err_count+=1
        if err_count>0:
            print('Could not load {} sequences'.format(err_count))
    return (seqs,np.array(xpos),np.array(ypos))

### CIFS ###

def retrieve_cif_names(run_path,cycle_nums):
    """
    Retrieve a list of cycle names
    
    Input:
    ---------------------
    run_path: the top directory
    cycle_nums: which cycle numbers to retrieve
    
    Outputs:
    ---------------------
    cycle_list: array of all the cycle names
    """
    all_cifs=[os.path.basename(filename) for filename in glob.glob(os.path.join(run_path,'cifs','*'))]
    cycle_cif={}
    for cif_file in all_cifs:
        cycle_string=cif_file.split('.')[0]
        cycle_cif[cycle_string]=cif_file
    cycle_list=[cycle_cif['C{}'.format(cycle_n)] for cycle_n in cycle_nums]
    return cycle_list

def gen_cif_filename(run_path,cycle_name,tile_num):
    """
    Generates the filename for the fastq file
    
    Input:
    ---------------------
    run_path: file path to run data
    cycle_name: name of cycle
    tile_num: the tile to use
    
    Outputs:
    ---------------------
    cif_file: path and file name for cif
    
    """
    cif_tile_file_format='s_1_{}.cif'.format(tile_num)
    cif_file=os.path.join(run_path,'cifs',cycle_name,cif_tile_file_format)
    
    return cif_file


def parse_cif(cif_file,n_channels=4):
    """
    Parse cif file into intensities by channel: A,C,T,G (or custom)
    
    Input:
    ---------------------
    cif_file: String with file and path for cif file
    n_channels: Number of channels in the file
    
    Outputs:
    ---------------------
    intensities: Intensities for the cif file of dim (n_channels,n_clusters)
    
    """
    with open(cif_file, 'rb') as file:
        # Extract all info
        cifs_identifier=np.fromfile(file,dtype=np.ubyte,count=3) # identifier "CIF"
        cifs_version=np.fromfile(file,dtype=np.int8,count=1) # version
        cifs_datatype=np.fromfile(file,dtype=np.int8,count=1) # data type, # of bytes
        cifs_firstcycle=np.fromfile(file,dtype=np.int16,count=1) # first cycle in file
        cifs_numcycles=np.fromfile(file,dtype=np.int16,count=1) # number of cycles
        cifs_numclusters=np.fromfile(file,dtype=np.int32,count=1) # number of clusters
        cifs_clusterintensities=np.fromfile(file,dtype='<H') # cluster intensities
    #     cifs_clusterintensities=np.fromfile(file,dtype=np.int16) # cluster intensities

    n_clusters=cifs_numclusters[0]
    # By channel: A,C,G,T
    intensities=cifs_clusterintensities[0:n_clusters*n_channels].reshape([n_channels,n_clusters])
    
    return intensities

### LINK FASTQ x/y to LOCS x/y


def round_locs_converted(locs_pts_converted,round_func1=np.round,round_func2=np.round):
    """
    Round converted locs vis the input functions
    
    Input:
    ---------------------
    locs_pts_converted: locs points to be rounded
    round_func1, round_func2: rounding functions
    
    Outputs:
    ---------------------
    locs_pts_rounded: points scaled, shifted and (now) rounded
    """
    locs_pts_rounded=np.zeros_like(locs_pts_converted)
    locs_pts_rounded[0]=round_func1(locs_pts_converted[0])
    locs_pts_rounded[1]=round_func2(locs_pts_converted[1])
    return locs_pts_rounded.astype(int)

def gen_locs_candidates(locs_pts,scale=10,shift=1000):
    """
    Convert locs to fastq range via shifting, scaling, and rounding
    
    Input:
    ---------------------
    locs_pts: Points to be changed
    round_func1, round_func2: rounding functions
    scale, shift: default values
    
    Outputs:
    ---------------------
    locs_pts_converted: points scaled, shifted and rounded
    """
    # Scale and shift
    locs_pts_converted=scale*locs_pts+shift
    
    # Prepare to store rounded points
    n_dim,n_pts=locs_pts_converted.shape
    n_round_comb=5
    locs_converted_rounded=np.zeros((n_round_comb,n_dim,n_pts),dtype=int)
    
    # Compute and store each combination of rounding
    func1,func2=np.round,np.round
    locs_converted_rounded[0,:,:]=round_locs_converted(locs_pts_converted,round_func1=func1,round_func2=func2)
    for i in range(2):
        func1=np.floor if i==0 else np.ceil
        for j in range(2):
            func2=np.floor if j==0 else np.ceil
            locs_converted_rounded[i*2+j+1,:,:]=round_locs_converted(locs_pts_converted,round_func1=func1,round_func2=func2)
    
    return locs_converted_rounded

def convert_1d(corr_max,pts):
    """
    Convert set of input pts to 1D form
    
    Input:
    ---------------------
    corr_max: max values for each dimensions (n_dim=2)
    pts: points to be converted (n_dim=2,n_pts)
    
    Outputs:
    ---------------------
    pts_1d: points in 1d form
    pts_1d_sort_ind: indices of sorted (argsort)
    input_1d: sorted 1d points
    """
    pts_1d=np.array([i+corr_max[0]*j for (i,j) in pts.T])
    pts_1d_sort_ind=np.argsort(pts_1d)
    input_1d=pts_1d[pts_1d_sort_ind]
    
    return (pts_1d,pts_1d_sort_ind,input_1d)

def intersect_locs_fastq(corr_max,locs_1d_data,fastq_1d_data):
    """
    Find intersection between two 1D sets
    
    Input:
    ---------------------
    corr_max: max values for each dimensions (n_dim=2)
    locs_1d_data: Contains 1D pts, indices, and sorted pts for locs
    fastq_1d_data: Contains 1D pts, indices, and sorted pts for fastq
    
    Outputs:
    ---------------------
    locs_linked_ind: linked indices for original locs of points found
    fastq_linked_ind: linked indices for original fastq of points found
    fastq_pts_found: bool array of points found
    """
    # Extract 1D data
    fastq_pts_1d,fastq_pts_1d_sort_ind,fastq_1d=fastq_1d_data
    locs_pts_1d,locs_pts_1d_sort_ind,locs_1d=locs_1d_data
    
    # Find loc points overlap
    pts_overlap=np.isin(locs_1d,fastq_1d)
    pts_overlap_ind=np.where(pts_overlap)[0] # Use 0 or returns in a list otherwise

    # Reverse find fastq points
    fastq_pts_overlap=np.isin(fastq_1d,locs_1d[pts_overlap_ind])
    fastq_pts_overlap_ind=np.where(fastq_pts_overlap)[0]
    
    # Use found values to retrieve indices on original arrays
    locs_linked_ind=locs_pts_1d_sort_ind[pts_overlap]
    fastq_linked_ind=fastq_pts_1d_sort_ind[fastq_pts_overlap]

    # Mark points that have been found
    fastq_pts_found=np.zeros_like(fastq_1d,dtype=bool)
    fastq_pts_found[fastq_linked_ind]=True
    
    return (locs_linked_ind,fastq_linked_ind,fastq_pts_found)


def link_fastq_to_ind(n_fastq_seqs,locs_converted_rounded,corr_max,fastq_1d_data):
    """
    Link fastq index to loc/cif index
    
    Input:
    ---------------------
    n_fastq_seqs: number of total fastq sequences in tile
    locs_converted_rounded: rounded locs data (combinations of round, floor, ceil)
    corr_max: max values for each dimensions (n_dim=2)
    fastq_1d_data: Contains 1D pts, indices, and sorted pts for fastq
    
    Outputs:
    ---------------------
    fastq_seq_ind: linked fastq to loc/cif index
    fastq_pts_found_total: bool array of which sequences were found
    """
    # Set up total arrays for found and seq-index
    fastq_seq_ind=-1*np.ones((n_fastq_seqs),dtype=int)
    fastq_pts_found_total=np.zeros((n_fastq_seqs),dtype=bool)

    # Search through each rouncding combination
    n_round_comb=locs_converted_rounded.shape[0]
    for i in range(n_round_comb):
        locs_pts_converted_i=np.squeeze(locs_converted_rounded[i])
        locs_1d_data_i=convert_1d(corr_max,locs_pts_converted_i)
        (locs_linked_ind_i,fastq_linked_ind_i,fastq_pts_found_i)=intersect_locs_fastq(
            corr_max,locs_1d_data_i,fastq_1d_data)

        # Link sequences
        fastq_seq_ind[fastq_linked_ind_i]=locs_linked_ind_i

        # Update total found
        fastq_pts_found_total=np.any(np.vstack((fastq_pts_found_total,fastq_pts_found_i)),axis=0)

    #     print('[{}] Total linked: {} of {}'.format(i,np.sum(fastq_pts_found_total),len(fastq_pts_found_total)))

    return (fastq_seq_ind,fastq_pts_found_total)

### SEQUENCE FILTERING ###

def find_seq_comp(seq):
    """
    Converts sequence to complement
    
    Input:
    ---------------------
    seq: sequence to convert
    
    Outputs:
    ---------------------
    seq_comp: Complement of sequence
    """
    base_comp={'A':'T','T':'A','C':'G','G':'C'}
    seq_comp=''
    for base in seq:
        seq_comp+=base_comp[base]
    return seq_comp

def find_seq_rev_comp(seq):
    """
    Converts sequence to reverse-complement
    
    Input:
    ---------------------
    seq: sequence to convert
    
    Outputs:
    ---------------------
    seq_rev_comp: Reverse complement of sequence
    """
    return find_seq_comp(seq[::-1])

def gen_seq_regex(input_info,seq_format='single'):
    """
    Generates the seq_regex format
    
    Input:
    ---------------------
    input_info: contains the format to be used
    seq_format: the format of the regex ('single' for single sequence, 'primers' for primer flanked)
    
    Outputs:
    ---------------------
    seq_regex: regex of the defined form
    """
    seq_regex=''
    
    # Operate based on inputs
    if seq_format=='primers':
        fp,rp=input_info
        rp_rev_comp=find_seq_rev_comp(rp)
        seq_regex='{}([ACTG]*){}'.format(fp,rp_rev_comp)
    else:
        seq_regex='({})'.format(input_info)
        
    return seq_regex

def filter_seqs(seqs,seq_regex,seq_lengths=[],filter_by_len=True):
    """
    Filters sequences to find the random regions
    Criteria defined by regex ('<fp>([ACTG]*)<rp_rev_comp>')
    Length defined by array
    
    Input:
    ---------------------
    seqs: list of sequences to filter
    seq_regex: regex with forward primer and reverse primer (rev comp)
    seq_length: range of sequences with appropriate lengths
    filter_by_len: bool of whether to filter by length (default: True)
    
    Outputs:
    ---------------------
    seqs_filtered: random regions passing filter
    inds_filtered: indices from original array
    """
    # Import regex search
    import re
    
    # Prepare empty lists
    seqs_filtered,inds_filtered=[],[]
    len_filtered=None
    if filter_by_len:
        len_filtered=[]
    
    # Loop through and filter
    for i,seq in enumerate(seqs):
        seq_found=re.findall(seq_regex,seq)
        if not not seq_found:
            random_region=seq_found[0]
            seqs_filtered.append(random_region)
            inds_filtered.append(i)
            if filter_by_len:
                len_filtered.append(len(random_region))
            
    # Convert to arrays
    seqs_filtered=np.array(seqs_filtered)
    inds_filtered=np.array(inds_filtered)
    
    # Filter by length
    lengths_pass=np.ones_like(seqs_filtered,dtype=bool)
    if filter_by_len:
        seq_len_min,seq_len_max=seq_lengths
        lengths_pass=np.all(np.vstack((np.array(len_filtered)>seq_len_min,np.array(len_filtered)<seq_len_max)),axis=0)
    
    # Return
    return (seqs_filtered[lengths_pass],inds_filtered[lengths_pass])

### DATA CONNECTION AND EXPORTING ###

def write_seq_intensity_complete(tile_num,seqs,x,y,intensities,cycle_list,run_path,fastq_name,is_filtered=False):
    """
    Write sequences and intensities by channel: A,C,T,G
    
    Input:
    ---------------------
    tile_num: tile to write out
    seqs: list of sequences
    intensities: list of intensities
    cycle_list: cycles to write
    run_path: run path of the imager data
    fastq_name: which fastq to write
    
    Outputs:
    ---------------------
    CSV files written out as 'sequence, tile, cycle1,...,cyclen' (header included)
    
    """
    # Make holder directory if not present already
    intensities_dir=os.path.join(run_path,'int_files')
    if not os.path.exists(intensities_dir):
        os.makedirs(intensities_dir)
    
    # Write using each cycle
    n_cycles=len(cycle_list)
    channels=['A','C','T','G']
    for i,channel in enumerate(channels):
        file_name=os.path.join(intensities_dir,fastq_name+'_{}.csv'.format(channel))
        if is_filtered:
            file_name=os.path.join(intensities_dir,fastq_name+'_filt_{}.csv'.format(channel))
        
        # Create new file with header if does not exist
        if not os.path.exists(file_name):
            with open(file_name,'w') as file:
                file.write('seq,tile,x,y,{}\n'.format(','.join(cycle_list)))
        
        # Write values
        with open(file_name,'a') as file:
            for j,seq in enumerate(seqs):
                file.write(seq+',{},{},{},{}\n'.format(tile_num,x[j],y[j],','.join(map(str,intensities[:,i,j]))))
                
def write_fastq_intensities(cycle_list,tile_list,run_path,fastq_list,regex_input=None,filter_output=False):
    """
    Write sequences and intensities by fastq
    
    Input:
    ---------------------
    cycle_list: cycles to write
    tile_list: tile to process
    run_path: run path of the imager data
    fastq_list: list of fastq names to process
    regex_input: filtering input if needed
    filter_output: whether to filter
    
    Outputs:
    ---------------------
    CSV files written out as 'sequence, tile, cycle1,...,cyclen' (header included)
    
    """
    n_cycles=len(cycle_list)
    
    # Unpack filter params
    regex_formats,regex_seqs,seq_lengths=None,None,None
    if filter_output:
        regex_formats,regex_seqs,seq_lengths=regex_input
    
    for tile_num in tile_list:
        print('Working on tile: {}'.format(tile_num))
        # clock
        t_start=time.time()

        # Import LOCS
        locs_file=gen_locs_filename(run_path,tile_num)
        locs_xpos,locs_ypos=parse_locs_file(locs_file)
        locs_pts=np.vstack((locs_xpos,locs_ypos))

        # Prepare to link
        locs_converted_rounded=gen_locs_candidates(locs_pts)
        (n_round_comb,n_dim,n_pts)=locs_converted_rounded.shape
        corr_max=np.amax(locs_converted_rounded,axis=(0,2))
        print(' - locs imported')

        # Retrieve intensities
        int_all=np.zeros((n_cycles,4,n_pts),dtype=int)
        for i,cycle_name in enumerate(cycle_list):
            cif_file=gen_cif_filename(run_path,cycle_name,tile_num)
            intensities=parse_cif(cif_file)
            int_all[i,:,:]=intensities
    #         print('Finished cycle {}'.format(cycle_name))
        print(' - intensities imported')

        # Loop through fastq
        for k,fastq_name in enumerate(fastq_list):
            tile_file=gen_fastq_filename(run_path,fastq_name,tile_num)
            (seqs,fastq_xpos,fastq_ypos)=parse_tile_fastq(tile_file)
            fastq_pts=np.vstack((fastq_xpos,fastq_ypos))

            fastq_1d_data=convert_1d(corr_max,fastq_pts)

            # Find linked points
            n_fastq_seqs=fastq_pts.shape[1]
            (fastq_seq_ind,fastq_pts_found_total)=link_fastq_to_ind(n_fastq_seqs,locs_converted_rounded,corr_max,fastq_1d_data)
            
            # Filter sequences
            if filter_output:
                seq_regex=gen_seq_regex(regex_seqs[k],seq_format=regex_formats[k])
                filter_seq_len=False if regex_formats[k]=='single' else True
                seqs_filtered,inds_filtered=filter_seqs(seqs,seq_regex,seq_lengths[k],filter_by_len=filter_seq_len)
                
                # Check to see if sequences are still left
                if len(inds_filtered)==0:
                    continue
                write_seq_intensity_complete(tile_num,seqs_filtered,fastq_xpos[inds_filtered],fastq_ypos[inds_filtered],
                                             int_all[:,:,fastq_seq_ind[inds_filtered]],cycle_list,run_path,
                                             fastq_name,is_filtered=True)
            else:
                # Write out
                write_seq_intensity_complete(tile_num,seqs,fastq_xpos,fastq_ypos,
                                             int_all[:,:,fastq_seq_ind],cycle_list,
                                             run_path,fastq_name,is_filtered=False)
                
        print(' - files written')
        t_end=time.time()
        print(' - time: {}'.format(t_end-t_start))