import csv
import os
import re
import time
from typing import Literal

import pandas as pd
from Bio import Entrez

from ..step_build import stage_builder

RANKS = [
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
]


def get_empty_rank():
    return dict.fromkeys(RANKS, "")


def get_lineage_from_entrez(taxid, logger):
    try:
        handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
        records = Entrez.read(handle)
        handle.close()

        if not records or len(records) == 0:
            return None

        lineage_ex = records[0].get("LineageEx", [])
        rank_info = {}

        for tax_level in lineage_ex:
            rank = tax_level.get("Rank")
            name = tax_level.get("ScientificName")
            if rank and name and (RANKS is None or rank in RANKS):
                rank_info[rank] = name

        current_rank = records[0].get("Rank")
        current_name = records[0].get("ScientificName")
        if (
            current_rank
            and current_name
            and (RANKS is None or current_rank in RANKS)
        ):
            rank_info[current_rank] = current_name

        return rank_info
    except Exception as e:
        logger.warning(f"Error retrieving taxonomy for taxid {taxid}: {e}")
        return None


class AddLineageStage(stage_builder.StageBuilder):
    def __init__(
        self,
        config,
        heading=os.path.basename(__file__),
        in_dir="",
        out_dir="",
        in_suffix=".csv",
        out_suffix=".csv",
        lineage_db: str | Literal["nucleotide"] = "nucleotide",
        entrez_email: str = None,
    ):
        super().__init__(
            heading=heading, config=config, in_dir=in_dir, out_dir=out_dir
        )
        self.in_suffix = in_suffix
        self.out_suffix = out_suffix
        self.lineage_db = lineage_db
        self.parse_lineage_db(entrez_email)

    def parse_genus2otherlv(self, lineage_path: str):
        self.genus2otherlv = {}
        with open(lineage_path) as in_handle:
            reader = csv.DictReader(in_handle)
            for row in reader:
                genus_name = row["genus_name"]
                otherlv = [
                    row["kingdom_name"],
                    row["phylum_name"],
                    row["class_name"],
                    row["order_name"],
                    row["family_name"],
                ]
                self.genus2otherlv[genus_name] = otherlv

    def parse_lineage_db(self, entrez_email):
        if self.lineage_db is None:
            pass
        elif os.path.isfile(self.lineage_db):
            self.parse_genus2otherlv(self.lineage_db)
            self.get_lineage_func = self.get_lineage_from_custom_db
        elif self.lineage_db == "nucleotide":
            self.get_lineage_func = self.get_lineage_from_entrez
            if entrez_email is None:
                self.entrez_email = "your.email@example.com"
                self.config.logger.warning(
                    "No email provided for Entrez to remotely query lineages. Set to 'your.email@example.com'. This email is used by NCBI to contact you in case of excessive usage or issues with your queries."
                )
            else:
                self.entrez_email = entrez_email
        else:
            raise ValueError(
                "lineage_db must be a valid file path or 'nucleotide'."
            )

    def setup(self, prefix):
        self.infile = os.path.join(self.in_dir, f"{prefix}{self.in_suffix}")
        self.outfile = os.path.join(self.out_dir, f"{prefix}{self.out_suffix}")
        self.check_infile()
        if not self.config.dry and hasattr(self, "get_lineage_func"):
            super().add_stage_function(
                "Add taxonomy to BLAST result",
                self.get_lineage_func,
            )

    def get_lineage_from_custom_db(self):
        if not hasattr(self, "genus2otherlv"):
            self.config.logger.warning(
                f"lineage_db is not set: {self.lineage_db}"
            )
            return False

        blast_result = pd.read_csv(self.infile)

        sseqids = blast_result["sseqid"].unique()

        lineage_data = []
        for sseqid in sseqids:
            lineage_record = {"sseqid": sseqid}
            species = sseqid.split("|")[-1]
            genus = species.split("_")[0]
            if genus in self.genus2otherlv.keys():
                taxonomy_levels = self.genus2otherlv[genus] + [genus, species]
            elif genus.endswith("idae") or genus.endswith("inae"):
                for _, otherlv in self.genus2otherlv.items():
                    if otherlv[-1] == genus:
                        taxonomy_levels = otherlv + ["", species]
                        break
            elif genus.endswith("iformes") or genus.endswith("oidea"):
                for _, otherlv in self.genus2otherlv.items():
                    if otherlv[-2] == genus:
                        taxonomy_levels = otherlv[:-1] + ["", "", species]
                        break
            else:
                self.config.logger.error(
                    f"genus {genus} not found in: {self.lineage_db}"
                )
                taxonomy_levels = [""] * 5 + [genus, species]

            lineage_record.update(
                dict(zip(RANKS, taxonomy_levels, strict=False))
            )
            lineage_data.append(lineage_record)

        lineage_df = pd.DataFrame(lineage_data)

        merged_result = blast_result.merge(
            lineage_df, left_on="sseqid", right_on="sseqid", how="left"
        )

        merged_result.to_csv(self.outfile, index=False)

        return True

    def get_lineage_from_entrez(self):
        blast_result = pd.read_csv(self.infile)

        Entrez.email = self.entrez_email

        sseqids = blast_result["sseqid"].unique()

        p_acc_id = re.compile(pattern=r"[A-Z]+_?\d+(?:\.\d+)?")

        lineage_data = []
        for sseqid in sseqids:
            lineage_record = {"sseqid": sseqid}

            m_acc_id = re.search(p_acc_id, sseqid)
            if m_acc_id:
                acc_id = m_acc_id.group()
            else:
                lineage_record.update(get_empty_rank())
                lineage_data.append(lineage_record)
                self.config.logger.warning(
                    f"No accession ID found for {sseqid}"
                )
                continue

            try:
                seqio = Entrez.efetch(
                    db=self.lineage_db, id=acc_id, retmode="xml"
                )
                seqio_read = Entrez.read(seqio)
                seqio.close()

                if not seqio_read or len(seqio_read) == 0:
                    lineage_record.update(get_empty_rank())
                    lineage_data.append(lineage_record)
                    self.config.logger.warning(
                        f"No information found for {acc_id}"
                    )
                    continue

                feature_table = seqio_read[0].get("GBSeq_feature-table", [])
                taxid = None

                for feature in feature_table:
                    if feature.get("GBFeature_key") == "source":
                        for qual in feature.get("GBFeature_quals", []):
                            if qual.get("GBQualifier_name") == "db_xref":
                                db_xref = qual.get("GBQualifier_value", "")
                                if db_xref.startswith("taxon:"):
                                    taxid = db_xref.split(":")[1]
                                    break

                if taxid:
                    rank_info = get_lineage_from_entrez(
                        taxid, self.config.logger
                    )
                    if rank_info:
                        for rank in RANKS:
                            lineage_record[rank] = rank_info.get(rank, "")
                    else:
                        lineage_record.update(get_empty_rank())
                        self.config.logger.warning(
                            f"No information found for {acc_id} (taxid: {taxid})"
                        )
                else:
                    lineage_record.update(get_empty_rank())
                    self.config.logger.warning(
                        f"No taxon ID found for {acc_id}"
                    )

                lineage_data.append(lineage_record)
                time.sleep(0.5)

            except Exception as e:
                lineage_record.update(get_empty_rank())
                lineage_data.append(lineage_record)
                self.config.logger.warning(f"Error processing {acc_id}: {e}")

        lineage_df = pd.DataFrame(lineage_data)

        merged_result = blast_result.merge(
            lineage_df, left_on="sseqid", right_on="sseqid", how="left"
        )

        merged_result.to_csv(self.outfile, index=False)

        return True

    def run(self):
        super().run()
        return all(self.output)


def addlineage_demo(
    config,
    prefix,
    blast_dir="",
    save_dir="",
    lineage_path="",
):
    stage = AddLineageStage(
        config,
        in_dir=blast_dir,
        out_dir=save_dir,
        lineage_db=lineage_path,
    )
    stage.setup(prefix)
    is_complete = stage.run()
    return is_complete
