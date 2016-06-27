# Script used to extract text from a collecion of tweets. 
# Tokenized and stemmed tweets are collected in a .gz file, one tweet per row. 
# We also use a binary file to store tweets ids, 8 bytes each.

import sys
import codecs
import struct
import json,re,os,sys
import gzip
from nltk.tokenize import word_tokenize
from unidecode import unidecode
from nltk.stem.snowball import EnglishStemmer
from nltk.corpus import stopwords
STOPWORDS = set(stopwords.words('english'))

emoticons_str = r"""
    (?:
        [:=;] # Eyes
        [oO\-]? # Nose (optional)
        [D\)\]\(\]/\\OpP] # Mouth
    )"""

regex_str = [
    emoticons_str,
    r'<[^>]+>', # HTML tags
    r'(?:@[\w_]+)', # @-mentions
    r"(?:\#+[\w_]+[\w\'_\-]*[\w_]+)", # hash-tags
    r'http[s]?://(?:[a-z]|[0-9]|[$-_@.&+]|[!*\(\),]|(?:%[0-9a-f][0-9a-f]))+', # URLs
    r'(?:(?:\d+,?)+(?:\.?\d+)?)', # numbers
    r"(?:[a-z][a-z'\-_]+[a-z])", # words with - and '
    r'(?:[\w_]+)', # other words
    r'(?:\S)' # anything else
]

tokens_re = re.compile(r'('+'|'.join(regex_str)+')', re.VERBOSE | re.IGNORECASE)
emoticon_re = re.compile(r'^'+emoticons_str+'$', re.VERBOSE | re.IGNORECASE)

def filter_words(words):
       return filter(lambda w: len(w)>=3 and w not in STOPWORDS, words)

def stem_words(words):
        s = EnglishStemmer()
        return map(s.stem, words)

def normalize_words(words):
        return map(lambda w: unidecode(w.lower()), words)

print "Usage:", sys.argv[0], "Directory OutputFile" 

def tokenize(s):
    return tokens_re.findall(s)
 
def preprocess(s, lowercase=False):
    tokens = tokenize(s)
    if lowercase:
        tokens = [token if emoticon_re.search(token) else token.lower() for token in tokens]
    filtered_t = filter_words(normalize_words(tokens))
    stemmed = stem_words(filtered_t)
    return stemmed

fids = open(sys.argv[2]+".ids", "wb")


with gzip.GzipFile(sys.argv[2]+".tweets.gz", 'w') as myzip:
  for root, dirs, files in os.walk(sys.argv[1]):
     for file in files:
        print file
        if 'json.gz' in file:
           print "Processing file: ", file
           zipfr = gzip.open("{}".format(os.path.join(root, file)), mode='rb')
           for line in zipfr:
              tweet = json.loads(line)
              tokens = preprocess(tweet['text'],True)
              if len(tokens) < 4:
                 continue
              text = " ".join(tokens)
              # print text       
              myzip.write(text + "\n")
              fids.write(struct.pack('<Q', long(tweet['id_str'])))

fids.close()
