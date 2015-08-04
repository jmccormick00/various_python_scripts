'''
Created on Feb 11, 2014

@author: jmccormick
'''

import pyodbc
#import pandas as pd
import pandas.io.sql as pdsql

if __name__ == '__main__':
    con = pyodbc.connect('DRIVER={SQL Server};SERVER=SERVERNAME;DATABASE=DATABASENAME;UID=auodbc;PWD=PASSWORD')
    
    sql = "select * from BLAH"
    
    masterDF = pdsql.frame_query(sql, con)
    con.close()
