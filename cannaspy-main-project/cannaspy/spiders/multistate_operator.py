import json
import scrapy
from scrapy.http import JsonRequest


class MultiStateOperatorSpider(scrapy.Spider):
    name = "multistate"
    allowed_domains = ["cannaspyglass.com"]
    headers = {
        "authority": "api.cannaspyglass.com",
        "accept": "*/*",
        "accept-language": "en-US,en;q=0.9",
        "origin": "https://portal-prd.cannaspyglass.com",
        "referer": "https://portal-prd.cannaspyglass.com/",
    }

    def start_requests(self):
        url = "https://api.cannaspyglass.com/csgExecStoredProc/?pname=csgGetMultiStateOperator&qstring=nathan@vernebio.com,name"
        yield JsonRequest(url=url, headers=self.headers, callback=self.parse_name)

    def parse_name(self, response):
        data = json.loads(response.text).get("data", [])
        if not data:
            return

        for business in data:
            holding_company_name = business["Business Name"]
            url = f"https://api.cannaspyglass.com/csgExecStoredProc/?pname=csgGetMultiStateOperatorDetails&qstring=nathan@vernebio.com%7Cname%7C{holding_company_name}"
            yield JsonRequest(
                url=url,
                headers=self.headers,
                meta={
                    "company_name": holding_company_name,
                    "num_licences": business["No Of Licenses"],
                    "num_state": business["No Of States"],
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
            item["No Of States"] = response.meta["num_state"]
            yield item
