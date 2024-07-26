import json
import scrapy
from scrapy.http import JsonRequest


class AssociationSpider(scrapy.Spider):
    name = "association"
    allowed_domains = ["cannaspyglass.com"]

    def start_requests(self):
        headers = {
            "authority": "api.cannaspyglass.com",
            "accept": "*/*",
            "accept-language": "en-US,en;q=0.9",
            "origin": "https://portal-prd.cannaspyglass.com",
            "referer": "https://portal-prd.cannaspyglass.com/",
        }
        url = "https://api.cannaspyglass.com/csgExecStoredProc/?pname=csgGetDataByAuthstates&qstring=nathan@vernebio.com%7Cvw_cannabisassociation%7Cstatename"
        yield JsonRequest(url=url, headers=headers, callback=self.parse_data)

    def parse_data(self, response):
        data = json.loads(response.text).get("data", [])
        if not data:
            return

        for item in data:
            yield item
