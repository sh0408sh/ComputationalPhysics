import requests

url = "https://yahoo-finance15.p.rapidapi.com/api/yahoo/ne/news"

headers = {
	"X-RapidAPI-Key": "06ac5ca0e8msh54d5b1db003ae49p10733bjsn41a386c3dcc4",
	"X-RapidAPI-Host": "yahoo-finance15.p.rapidapi.com"
}

response = requests.request("GET", url, headers=headers)

print(response.text)

