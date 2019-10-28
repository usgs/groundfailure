import gfail.pdl as pdl


def gf_transfer(event_dir, version=1, pdl_config=None, dry_run=False,
                status='UPDATE'):
    """
    Transfer ground failure results to dev server.

    Args:
        event_dir (str): Path to event directory.
        version (int): Version number to assign to ground-failure run.
        pdl_config (str): Path to PDL config file.
        dry_run (bool): True suppresses transfer but data is assembled for
            the transfer.
        status (str): Status of ground-failure product being sent to comcat.
            Default is "UPDATE" but can also be "WARNING" so that the product
            page displays the warning banner.

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

        log = pdl.transfer(
            event_dir, version, pdl_config,
            dryrun=dry_run, status=status)

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
