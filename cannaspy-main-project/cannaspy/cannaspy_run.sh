
output_folder="../results"

#scrapy crawl vendor -o ${output_folder}/vendor_data.csv

#scrapy crawl operator -o ${output_folder}/operator.csv

#scrapy crawl association -o ${output_folder}/association.csv

scrapy crawl multistate -o ${output_folder}/multistate_operator.csv

#scrapy crawl vertical -o ${output_folder}/vertical_data.csv
