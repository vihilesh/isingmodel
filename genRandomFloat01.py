import numpy as np
import base64
import json
from math import floor
import os
import requests





# randomHexBlock = "f87a9005dd5d1518d1fd6e539b0b14d5cc791b0bc977ce96ecd5f38bb14df431c4a64f57686758286039a2fc36d458ec3b87727e7839c8d31d86030e1d2014aa400714d1f7271b9b528f91251efac8cc219c6e867cab3594c97ac46840be9abffb38d5ec1998ef7db4fc6cde0d0fc5cd178a23e1e11438f12fb41a4b9ada5a2cd5bfe59b8955d0c25ae58b1dc65b03dbb4cd4e5e18072c84018efe56ae4276f0bd16dc85563eef0a8cb63217b68ed60e6f965f2d71dc4d0e0d0e1223be914510feb4e617f89faf45fad550b2aeb77e242a08fd994976210cd1001bd7b9a1dd5a520d3b2f98e414ffa778caad9f52647fddaba47cf126734f2dd9ba1cc26fa9e56926fee5c2962d8cd362bea50fe978278cc839850d97f43d57389f377010b3070710ecc17059e978ece5c1d97709ca535f4e3b41f3b993141c0486eff338256916c025e4a6c14c4594847c4ea0ed19c2a5051545905e3f3086f4a641333be0215ed44340d6342de8caf3b95f607e6158550f57bec8d235d775359167f727f14a3d7bf845a076282a1adf37e240009b757806f77894b17d9cd9da89c81b7ba10aa7047348635e6f87ee919ece5107d3cedfc42fcd74e76bd22edb02bff3d6cfbf06bf5972ad687fdaea3d88c93643046e8091a4c26fb0f3087879e6739bc9ca7c5b9bc2dd1e6b59da049cbcd6af7f236a64a3721c8317ca05773bb73545f92ce7ba6e595b626b154a657fe253e4fa360e255d05661a903e5f64775aa0e98198335b38e26bfcb92bd12fd26557836e2dc0be465cce2bab54e5eaaf45ad56338f86addf3568365a2f4034ef76c9c6a55fe2ded831ddf6c1dfe4ba8e53eb65b4fd0a45d21b09e806a28c0b5668655f573481347785259b315fb46e994457460b9abc9af58a2ff87f845e0705ecd809e39c3086169bbae9975943366d8ee1d61b9cbe4eebb7af2713a5b6bf2355a0a7f6f9207ed513dd2f42e121e86948374664cd0e73f8a7f74eeb7fbaf2cc15417c0481c9f9a92147b67cfc2035186dc857b97326ad75d4368cb980d4840514e26351ba382496dc7f2353d564f85d8ce02e5210875f14c4ad3de0ee98c4d66f1efcb3ab12392b98c0c03857e8ee43ccef8adf7387b41698ccd09a95201cb18f103b44c95d982d543156b8cb12af5c5c156b4d0496fbd4eaf60ded9c59232e8593330a729424d64965f27fe504577a644dd79b3983a6814b4dca8a43d57e771a29ccff6c02e8b99aaf8f3aa792b622a1275af9136fd87dd0833468346d0acd958e90fd7d42c5ac4fafa5ea1b754889f4851f666a61e1e68a616460f9d163a29a2a362a51de9de1643c0671fa77120148007f987343d3a4304136408d5593dadd9550d74c1ab39ae8f4104beb3ab0b61ba2d235088dd3f9275fb57eb0cdd94c23321ba827f71f29fb38c2da6e455f9bbfc36433ee7d"
# n = 6
# floatnumbers = []

def convert(number,points):
    decimal = pow(10,points) # power function

    return number / decimal

def convert(number):
    decimal = pow(2,64) # power function

    return number / decimal



headers = {
'Accept': '*/*',
'Accept-Encoding': 'gzip, deflate, br, zstd',
'Accept-Language': 'en-US,en;q=0.9',
'Connection': 'keep-alive',
'Content-Type': 'application/x-www-form-urlencoded; charset=UTF-8',
'Host':'qrng.anu.edu.au',
'Referer': 'https://qrng.anu.edu.au/random-block-hex/',
'Sec-Ch-Ua': '"Chromium";v="122", "Not(A:Brand";v="24", "Google Chrome";v="122"',
'Sec-Ch-Ua-Mobile': '?0',
'Sec-Ch-Ua-Platform': "Windows",
'Sec-Fetch-Dest': "empty",
'Sec-Fetch-Mode':"cors",
'Sec-Fetch-Site':'same-origin',
'User-Agent':'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/122.0.0.0 Safari/537.36',
'X-Requested-With':'XMLHttpRequest'
}





# The API endpoint to communicate with
url_get = "https://qrng.anu.edu.au/wp-content/plugins/colours-plugin/get_block_hex.php?_=1709402362392"

# A POST request to tthe API
get_response = requests.get(url_get, headers=headers, verify=False)
randomhexblock = get_response.text

hexfile = open("QRNGFloat64-8.txt", "a") 

for i in range(0, len(randomhexblock), 16): 
    substring = randomhexblock[i:i+16]
    temp = (int(substring,16))
    #templength = len(str(temp))
    #floatnum = convert(temp, templength)
    floatnum = convert(temp)
    hexfile.write(str(floatnum)+"\n")
   
hexfile.close


