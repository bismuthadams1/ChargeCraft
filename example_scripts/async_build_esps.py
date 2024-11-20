from chargecraft.storage.qcarchive_transfer import QCArchiveToLocalDB
from qcportal import PortalClient
from chargecraft.storage.storage import MoleculePropRecord, MoleculePropStore
from openff.recharge.grids import LatticeGridSettings
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import Pool, get_context
from tqdm import tqdm

def main():

    client = PortalClient("api.qcarchive.molssi.org")

    prop_store = MoleculePropStore("/mnt/storage/nobackup/nca121/QC_archive_50K_esp_gen/async_chargecraft/ESP_rebuilt.db")

    grid_settings = LatticeGridSettings(
    type="fcc", spacing=0.5, inner_vdw_scale=1.4, outer_vdw_scale=2.0
    )

    molecules = [record for record in client.query_records(dataset_id=347)]

    qcarchdb = QCArchiveToLocalDB(
        qc_archive=client,
        prop_data_store=prop_store,
        grid_settings=grid_settings,
    )


    with ProcessPoolExecutor(
        max_workers=2, mp_context=get_context("spawn")
    ) as pool:

        futures = [
            pool.submit(
                qcarchdb.process_item(item = molecule,
                                      exclude_keys = None,
                                      qm_esp = True,
                                      compute_properties = True,
                                      return_store = True
                                      )
            )
            for molecule in molecules
        ]
        #to avoid simultaneous writing to the db, wait for each calculation to finish then write
        for future in tqdm(as_completed(futures, timeout=300), total=len(futures)):
            esp_record = future.result()
            prop_store.store(esp_record)



if __name__ == "__main__":
    main() 