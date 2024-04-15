import sqlite3


class QueryExecutor:
    
    def __init__(self, db_conn, upstream, downstream):
        self.db_conn = db_conn #sqlite3.connect(db_path)
        self.upstream = upstream 
        self.downstream = downstream

    def query(self, sql, category, params=()):
        """
        Execute a query with the given SQL and parameters.

        :param sql: A string containing the SQL query to execute.
        :param params: A tuple of parameters to bind to the query.
        :return: The query results or None if an error occurred.
        """
        #####################
        ## TODO
        if category == 'id_convert':
            cur_params = (
                params[0],
                params[1],
            )
        #####################
        if category == 'posi2gene':
            
            cur_params = (
                params[0],
                params[1],
                int(params[2]) + self.upstream,
                int(params[2]) - self.downstream
            )

        #####################
           
        elif category in ['posi2exon','posi2feature','posi2cfeature1','posi2cfeature_sd','posi2cfeature_sa']:
            
            cur_params = (
                params[0],
                params[1],
                params[1]
            )

        elif category in ['test','ekid2geneid','getGenelimit','extract_geneset']:
            cur_params = (params)
        
        elif category in ['extract_tf_gene']:
            cur_params = ()
        
        #####################  
        try:
            cursor = self.db_conn.cursor()
            cursor.execute(sql, cur_params)
            results = cursor.fetchall()
            cursor.close()
            return results if results else None
        
        except sqlite3.Error as e:
            print(f"Database error: {e}")
            return None
        
        except ValueError as e:
            print(f"Error: {e}")
            return None

        
    def close(self):
        self.db_conn.close()