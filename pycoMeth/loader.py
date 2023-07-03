import logging

from meth5 import MetH5File

from typing import List, Tuple, Dict
from pycoMeth.CoordGen import Coord


def read_sample_ids_from_read_groups(h5_file_list, read_group_key, labels=None):
    rg_dict = {}
    for fn in h5_file_list:
        with MetH5File(fn, "r") as f:
            f_rg_dict = f.get_all_read_groups(read_group_key)
            for k, v in f_rg_dict.items():
                if k in rg_dict and rg_dict[k] != v:
                    raise pycoMethError("Read groups in meth5 files must have the same encoding")
            rg_dict.update(f_rg_dict)
    if labels is not None:
        return [k for k, v in rg_dict.items() if v in labels]
    else:
        return [k for k in rg_dict]


class MetH5Loader:
    def __init__(self, h5_read_groups_key: str, sample_id_list: List, h5_file_list: List):
        self.h5_read_groups_key = h5_read_groups_key
        self.sample_hf_files: Dict[str, MetH5File] = {}
        self.llr_threshold = 2.0  # TODO expose parameter
        
        if h5_read_groups_key is None:
            for sample_id, h5_file in zip(sample_id_list, h5_file_list):
                hf = MetH5File(h5_file, "r")
                self.sample_hf_files[sample_id] = hf
        else:
            hf = MetH5File(h5_file_list[0], "r")
            for sample_id in sample_id_list:
                self.sample_hf_files[sample_id] = hf
    
    def __del__(self):
        for hf in self.sample_hf_files.values():
            try:
                hf.close()
            except:
                pass
    
    def read_raw_llrs(self, interval: Coord) -> Tuple[List, List, List, List]:
        sample_llrs = {}
        sample_pos = {}
        sample_reads = {}
        for sample_id, hf in self.sample_hf_files.items():
            chrom_container = hf[interval.chr_name]

            if chrom_container is None:
                continue

            interval_container = chrom_container.get_values_in_range(interval.start, interval.end)
            
            if interval_container is None:
                continue
            llrs = interval_container.get_llrs()[:]
            pos = interval_container.get_ranges()[:, 0]
            read_names = interval_container.get_read_ids()[:]
            
            if self.h5_read_groups_key is not None:
                read_samples = interval_container.get_read_groups(self.h5_read_groups_key)
                mask = [rs == sample_id for rs in read_samples]
                llrs = llrs[mask]
                pos = pos[mask]
                read_names = read_names[mask]
            
            sample_llrs[sample_id] = llrs.tolist()
            sample_pos[sample_id] = pos.tolist()
            sample_reads[sample_id] = read_names.tolist()
        
        # Remove samples for which there is no data
        label_list = list([k for k in sample_llrs.keys() if len(sample_llrs[k]) > 0])
        raw_llr_list = [sample_llrs[k] for k in label_list]
        raw_pos_list = [sample_pos[k] for k in label_list]
        raw_read_list = [sample_reads[k] for k in label_list]
        return label_list, raw_llr_list, raw_pos_list, raw_read_list
    
    @staticmethod
    def interpret_sample_ids_from_arguments(sample_id_list, read_groups_key, h5_file_list):
        if sample_id_list is None:
            if read_groups_key is None:
                # If no sample id list is provided and no read group key is set
                # automatically define tests and maximal missing samples depending on number of files to compare
                return list(range(len(h5_file_list)))
            else:
                return read_sample_ids_from_read_groups(h5_file_list, read_groups_key)
        elif read_groups_key is not None:
            # H5 file stores groups as int
            return read_sample_ids_from_read_groups(h5_file_list, read_groups_key, labels=sample_id_list)
        else:
            return sample_id_list
