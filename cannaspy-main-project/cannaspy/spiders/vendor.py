import json
import scrapy
from scrapy.http import JsonRequest
from .utils import US_STATE_CODE


class VendorSpider(scrapy.Spider):
    name = "vendor"
    allowed_domains = ["cannaspyglass.com"]

    def start_requests(self):
        headers = {
            "authority": "api.cannaspyglass.com",
            "accept": "*/*",
            "accept-language": "en-US,en;q=0.9",
            "origin": "https://portal-prd.cannaspyglass.com",
            "referer": "https://portal-prd.cannaspyglass.com/",
        }
        for code in US_STATE_CODE:
            url = f"https://api.cannaspyglass.com/csgExecStoredProc/?pname=csgGetVendorDataByState&qstring=nathan@vernebio.com%7C{code}"
            yield JsonRequest(url=url, headers=headers, callback=self.parse_data)

    def parse_data(self, response):
        data = json.loads(response.text).get("data", [])
        if not data:
            return

        for vendor in data:
            yield vendor
