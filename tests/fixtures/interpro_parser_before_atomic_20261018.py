"""Frozen original InterPro full-map parser for exact serialization parity only.

Copied before the2026-10-18 atomic-output/failure-propagation repair. This is a
reference implementation, never an installation route or fixture acceptance.
"""
import gzip
import logging
import time
import xml.etree.ElementTree as ET
from collections import defaultdict
from typing import Dict, List, Tuple
logger = logging.getLogger(__name__)
GOTermInfo = Tuple[str, str, str]
InterProMapValue = Tuple[str, str, str, str, List[GOTermInfo], str]
FinalMapping = Dict[Tuple[str, str], List[InterProMapValue]]

# Original function SHA256: 27184ab500ea5ec999320b511b2682a967c3848c40057d9f8f6a6d83677db578
def parse_interpro_xml_for_mapping(xml_file_path: str) -> FinalMapping:
    """
    Parses a gzipped InterPro XML file (interpro.xml.gz) to extract mappings
    from member database signatures to InterPro entries, including GO terms
    and abstract text. (Single-process version).
    """
    logger.info(f"Starting single-process parsing of {xml_file_path} for FULL mapping...")
    start_time = time.time()
    member_to_interpro_map: FinalMapping = defaultdict(list)
    interpro_entry_count = 0

    try:
        with gzip.open(xml_file_path, 'rb') as f:
            context = ET.iterparse(f, events=('end',))
            _, root = next(context)

            for event, elem in context:
                if elem.tag == 'interpro':
                    interpro_entry_count += 1
                    ipr_id = None
                    try:
                        ipr_id = elem.get('id')
                        ipr_type = elem.get('type')
                        ipr_short_name = elem.get('short_name')
                        name_element = elem.find('name')
                        ipr_name = name_element.text.strip() if name_element is not None and name_element.text else ""

                        go_terms_list: List[GOTermInfo] = []
                        class_list = elem.find('class_list')
                        if class_list is not None:
                            for classification in class_list.findall('classification'):
                                if classification.get('class_type') == 'GO':
                                    go_id = classification.get('id')
                                    cat_elem = classification.find('category')
                                    desc_elem = classification.find('description')
                                    go_category = cat_elem.text.strip() if cat_elem is not None and cat_elem.text else ""
                                    go_description = desc_elem.text.strip() if desc_elem is not None and desc_elem.text else ""
                                    if go_id:
                                        go_terms_list.append((go_id, go_category, go_description))

                        abstract_text = ""
                        abstract_element = elem.find('abstract')
                        if abstract_element is not None:
                            abstract_text = " ".join(t for t in abstract_element.itertext() if t).strip()
                            abstract_text = ' '.join(abstract_text.split())

                        map_value: InterProMapValue = (
                            ipr_id or "UnknownIPR",
                            ipr_type or "UnknownType",
                            ipr_short_name or "UnknownShortName",
                            ipr_name,
                            go_terms_list,
                            abstract_text
                        )

                        member_list = elem.find('member_list')
                        if member_list is not None and ipr_id is not None:
                            for db_xref in member_list.findall('db_xref'):
                                db_name = db_xref.get('db')
                                db_key = db_xref.get('dbkey')
                                if db_name and db_key:
                                    map_key = (db_name.upper(), db_key)
                                    member_to_interpro_map[map_key].append(map_value)

                    except Exception as inner_e:
                         logger.exception(f"Error processing details for InterPro entry near count {interpro_entry_count} (ID: {ipr_id})")
                    finally:
                         elem.clear()

                    if interpro_entry_count % 10000 == 0:
                        elapsed = time.time() - start_time
                        logger.info(f"Processed {interpro_entry_count} entries (full map)... ({elapsed:.1f}s)")

        end_time = time.time()
        logger.info(f"\nFinished FULL mapping parsing.")
        logger.info(f"Processed {interpro_entry_count} InterPro entries.")
        logger.info(f"Created mapping dictionary with {len(member_to_interpro_map)} unique member signature keys.")
        logger.info(f"Total time: {end_time - start_time:.2f} seconds.")
        return dict(member_to_interpro_map)

    except FileNotFoundError:
        logger.error(f"Error: File not found at {xml_file_path}")
        return {}
    except ET.ParseError as e:
        logger.exception(f"Error: Failed to parse XML file")
        return {}
    except Exception as e:
        logger.exception(f"An unexpected error occurred during file processing")
        return {}
