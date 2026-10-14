"""
Range-header pagination, for IEDB

IEDB's query API (https://query-api.iedb.org/api/v1) is PostgREST, which pages
results through the HTTP `Range` request header rather than a query parameter,
and signals the end of data with a 416 status rather than an empty page. Kept
here rather than under Sources/ directly since only specificity sources page
this way so far.
"""

# Info
__author__ = 'Pramod Shinde'

# Imports
import logging

# Immvoke imports
from immvoke.Exceptions import HttpError

log = logging.getLogger(__name__)


def pageByRange(client, url, page_size=2000, params=None, headers=None):
    """
    Page a Range-header paginated endpoint.

    Arguments:
      client (HttpClient): the shared HTTP client.
      url (str): the endpoint's absolute URL.
      page_size (int): rows requested per page.
      params (dict): query parameters sent with every request, e.g. PostgREST
        column filters.
      headers (dict): additional headers sent with every request.

    Yields:
      list: one page of parsed JSON records. A page shorter than page_size,
      or a Content-Range total that has been reached, ends iteration.

    Raises:
      HttpError: if a page comes back with a status other than 200, 206 or 416.
    """
    base_headers = dict(headers or {})
    start = 0
    fetched = 0

    while True:
        end = start + page_size - 1
        page_headers = {**base_headers, 'Range-Unit': 'items',
                        'Range': '%d-%d' % (start, end)}
        response = client.get(url, params=params, headers=page_headers)

        # 416 past the last row is how this style of API says "no more data".
        if response.status_code == 416:
            return
        if response.status_code not in (200, 206):
            raise HttpError('GET %s returned HTTP %d'
                            % (url, response.status_code))

        batch = response.json()
        if batch:
            yield batch
        fetched += len(batch)

        # Content-Range: items 0-1999/42000 reports the total row count, when
        # the endpoint sends it, so pagination can stop as soon as every row
        # has arrived rather than always waiting for a short final page.
        content_range = response.headers.get('Content-Range', '')
        total_text = content_range.rsplit('/', 1)[-1] if '/' in content_range else ''
        if total_text.isdigit() and fetched >= int(total_text):
            return

        if len(batch) < page_size:
            return
        start += page_size
