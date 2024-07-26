import json
import scrapy
from scrapy.http import JsonRequest
from .utils import US_STATE_CODE


class VerticalSpider(scrapy.Spider):
    name = "vertical"
    allowed_domains = ["cannaspyglass.com"]
    headers = {
        "authority": "api.cannaspyglass.com",
        "accept": "*/*",
        "accept-language": "en-US,en;q=0.9",
        "origin": "https://portal-prd.cannaspyglass.com",
        "referer": "https://portal-prd.cannaspyglass.com/",
    }

    def start_requests(self):
        for code in US_STATE_CODE:
            url = f"https://api.cannaspyglass.com/csgExecStoredProc/?pname=csgGetVerticalIntegrationByState&qstring=nathan@vernebio.com,name,{code}"
            yield JsonRequest(
                url=url,
                headers=self.headers,
                meta={"state": code},
                callback=self.parse_name,
            )

    def parse_name(self, response):
        data = json.loads(response.text).get("data", [])
        if not data:
            return

        for business in data:
            holding_company_name = business["Business Name"]
            url = f"https://api.cannaspyglass.com/csgExecStoredProc/?pname=csgGetOperatorDataByState&qstring=nathan@vernebio.com%7C{holding_company_name}%7C{response.meta['state']}%7Cname"
            yield JsonRequest(
                url=url,
                headers=self.headers,
                meta={
                    "company_name": holding_company_name,
                    "num_licences": business["No Of Licenses"],
                    "num_verticals": business["No Of Verticals"],
                },
                callback=self.parse_data,
            )

    def parse_data(self, response):
        data = json.loads(response.text).get("data", [])
        if not data:
            return

        for item in data:
            item["Business Name"] = response.meta["company_name"]
            item["No Of Licenses"] = response.meta["num_licences"]
            item["No Of States"] = response.meta["num_verticals"]
            yield item
