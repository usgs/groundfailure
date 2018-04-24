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

            # Construct URL to index.html
            se_lines = log['se'].decode().split("\n")
            info_line = [l for l in se_lines if 'send complete Socket' in l][0]
            info = info_line.split(':')
            eventid = info[-2]
            source = info[-4]

            # URL for event's detail feed
            feed = ('https://dev01-earthquake.cr.usgs.gov/fdsnws/event/1/'
                    'query?eventid=[EVENTID]&format=geojson')
            feed = feed.replace('[EVENTID]', source+eventid)

            # Print to screen/logger
            print("Event detail JSON feed: \n%s " % feed)
        elif log['rc'] is True and dry_run is True:
            print("Dry run complete, no transfer attempted.")
            success = True
            feed = ''
        else:
            print('PDL transfer failed.')
            print(log['so'].decode())
            success = False
            feed = ''

        return success, feed
