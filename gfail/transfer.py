import gfail.pdl as pdl


def gf_transfer(event_dir, pdl_config=None, dry_run=False):
    """
    Transfer ground failure results to dev server.

    Args:
        event_dir (str): Path to event directory.
        pdl_config (str): Path to PDL config file.
        dry_run (bool): True suppresss transfer but data is assesmbled for
            the transfer.

    Returns:
    tuple:
        - Success: bool for whether transfer was successful.
        - feed: link to json feed for product.
    """

    print('Preparing directory to transfer %s...'
          % event_dir)
    pdl.prepare_pdl_directory(event_dir)

    if pdl_config is None:
        print('PDL directory prepared, no pdl_config '
              'provided so no files were sent')
        return (True, None, None)
    else:
        # Transfer
        if not dry_run:
            print('Transferring...')
        else:
            print('Constructing PDL command...')

        log = pdl.transfer(event_dir, pdl_config, dryrun=dry_run)

        if log['rc'] is True and dry_run is False:
            print('Successful PDL transfer.')
            success = True
        elif log['rc'] is True and dry_run is True:
            print("Dry run complete, no transfer attempted.")
            success = True
        else:
            print('PDL transfer failed.')
            print(log['so'].decode())
            success = False

        return success
