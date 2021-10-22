import json

import pandas as pd


class MMP:
    def __init__(
        self,
        mmp_core: str,
        mmp_smirks: str,
        mmp_molecule_names: list[str],
        mmp_molecule_smiles: list[str],
    ):
        self.__core = mmp_core.strip()
        self.__smirks = [mmp_molecule_smiles, mmp_smirks]
        self.__molecule_names = mmp_molecule_names
        self.__molecule_smiles = mmp_molecule_smiles
        self.__core_sizes = len([x for x in mmp_core if x.isalpha()])

    def get_core(self):
        return self.__core

    def get_variable_parts(self):
        return self.__smirks

    def get_molecule_names(self):
        return self.__molecule_names

    def get_molecule_smiles(self):
        return self.__molecule_smiles

    def get_core_size(self):
        return self.__core_sizes


class MMPSeries:
    def __init__(self, first_mmp: MMP):
        self.__core = first_mmp.get_core()
        self.__variable_parts = [first_mmp.get_variable_parts()]
        self.__molecule_names = list(first_mmp.get_molecule_names())
        self.__molecule_smiles = list(first_mmp.get_molecule_smiles())
        self.__core_size = first_mmp.get_core_size()

    def _add_variable_part(self, addition):
        self.__variable_parts.append(addition)

    def _add_molecule_name(self, addition):
        self.__molecule_names.append(addition)

    def _add_molecule_smiles(self, addition):
        self.__molecule_smiles.append(addition)

    def add_mmp(self, single_mmp: MMP):
        if self.__core == single_mmp.get_core():
            print(
                f"Adding mmp {single_mmp.get_molecule_names()} to series with core {single_mmp.get_core()} and our own {self.__core}"
            )
            self._add_variable_part(single_mmp.get_variable_parts())
            self._add_molecule_name(single_mmp.get_molecule_names()[0])
            self._add_molecule_name(single_mmp.get_molecule_names()[1])
            self._add_molecule_smiles(single_mmp.get_molecule_smiles()[0])
            self._add_molecule_smiles(single_mmp.get_molecule_smiles()[1])

    def check_mmp_compatibility(self, single_mmp: MMP):
        return self.__core == single_mmp.get_core()

    def get_core(self):
        return self.__core

    def get_variable_parts(self):
        return self.__variable_parts

    def get_molecule_smiles(self):
        return self.__molecule_smiles

    def get_molecule_names(self):
        return self.__molecule_names

    def get_mms_size(self):
        return len(set(self.__molecule_names))

    def get_molecule_sizes(self):
        return self.__core_size

    def save_mms_data(self, path_to_save):
        data_dict = {
            "core": self.__core,
            "smirks": self.__variable_parts,
            "smiles": self.__molecule_smiles,
            "names": self.__molecule_names,
            "size": self.get_mms_size(),
        }
        with open(path_to_save, "w") as f:
            json.dump(data_dict, f)


class MMPOutputHandler:
    def __init__(self, output_path, rdkit=True, mmpdb=False, minimum_core=None):
        self.path = output_path
        self.rdkit = rdkit
        self.mmpdb = mmpdb
        self.mmps = []
        self.series = []
        self.minimum_core = minimum_core

    def read_in_mmps_rdkit(self):
        with open(self.path.as_posix(), "r") as f:
            lines = f.readlines()

        for line in lines:
            if self.mmpdb:
                parts = line.split("\t")
            else:
                parts = line.split(",")
            new_mmp = MMP(
                mmp_core=parts[5],
                mmp_smirks=parts[4],
                mmp_molecule_names=[parts[2], parts[3]],
                mmp_molecule_smiles=[parts[0], parts[1]],
            )
            print(
                f"New MMP of molecules {new_mmp.get_molecule_names()} with core {new_mmp.get_core()} and variable_part {new_mmp.get_variable_parts()}"
            )
            if self.minimum_core == None or new_mmp.get_core_size() > self.minimum_core:
                self.mmps.append(new_mmp)

    def read_in_mmps_naomi(self):
        with open(self.path.as_posix(), "r") as f:
            lines = f.readlines()

        all_mmps = []
        for line in lines:
            if line.startswith("MatchedMolecularPair"):
                core = ""
                smirks = ""
                molecule_names = [
                    line.split(" ")[2].split("(")[1],
                    line.split(" ")[7].split("(")[1],
                ]
                molecule_smiles = [line.split(" ")[5], line.split(" ")[10]]
            elif line.startswith("with"):
                smirks = line.split(" ")[2]
            elif line.startswith("and"):
                core = line.split(" ")[2]
            elif line.startswith("Value"):
                new_mmp = MMP(
                    mmp_core=core,
                    mmp_smirks=smirks,
                    mmp_molecule_names=molecule_names,
                    mmp_molecule_smiles=molecule_smiles,
                )
                all_mmps.append(new_mmp)

    def convert_mmps_to_series(self):
        for mmp in self.mmps:
            found_series = False
            for series in self.series:
                if series.check_mmp_compatibility(mmp):
                    series.add_mmp(mmp)
                    found_series = True
            if not found_series:
                self.series.append(MMPSeries(mmp))
        self.series.sort(
            key=lambda x: (x.get_mms_size(), x.get_molecule_sizes()), reverse=True
        )

    def output_series(self):
        for series in self.series:
            print(
                series.get_mms_size(),
                series.get_molecule_names(),
                series.get_core(),
                series.get_molecule_sizes(),
            )

    def save_usable_series(self, df_list):
        if len(df_list) == 0:
            return
        folder_path = self.path.parent / self.path.stem
        if not folder_path.exists():
            folder_path.mkdir()
        for index, df in enumerate(df_list):
            final_path = folder_path / f"mmp_{index}.smi"
            data_path = folder_path / f"mmp_{index}_data.json"
            df[0].to_csv(final_path.as_posix(), index=None, header=None, sep=" ")
            df[1].save_mms_data(data_path)

    def run(self):
        if self.rdkit:
            self.read_in_mmps_rdkit()
        else:
            self.read_in_mmps_naomi()
        self.convert_mmps_to_series()
        usable_series = self.calc_usable_series()
        self.save_usable_series(usable_series)

    def calc_usable_series(self, cutoff=2):
        used_complexes = []
        usable_series = []
        for series in self.series:
            non_chembl_ligands = [
                x
                for x in series.get_molecule_names()
                if not (x.startswith("CHEMBL") or x in used_complexes)
            ]
            works = len(non_chembl_ligands) > 0 and series.get_mms_size() >= cutoff
            if works:
                usable_series.append(series)
                used_complexes.extend(non_chembl_ligands)
        usable_dataframes = []
        for series in usable_series:
            pre_df = list(
                zip(series.get_molecule_smiles(), series.get_molecule_names())
            )
            df = pd.DataFrame(pre_df)
            df.drop_duplicates(inplace=True)
            usable_dataframes.append([df, series])
        return usable_dataframes
