import json
import scrapy
from scrapy.http import JsonRequest
from .utils import US_STATE_CODE


class OperatorSpider(scrapy.Spider):
    name = "operator"
    allowed_domains = ["cannaspyglass.com"]

    def start_requests(self):
        headers = {
            "referer": "https://portal-prd.cannaspyglass.com/",
        }
        for code in US_STATE_CODE:
            url = f"https://portal-dev.cannaspyglass.com/prd/AppLoad/LicenseInfo/BusinessEntities_{code}.json"
            yield JsonRequest(url=url, headers=headers, callback=self.parse_data)

    def parse_data(self, response):
        data = json.loads(response.text)
        if not data:
            return

        for item in data:
            yield item
