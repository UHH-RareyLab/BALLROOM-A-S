import unittest

from src.utils.db_interaction import query_chembl, query_activitydb


class MyTestCase(unittest.TestCase):
    def test_query_chembl(self):
        simple_query = f"SELECT * from Version"
        data = query_chembl(simple_query)
        self.assertTrue(data.shape[0] > 0)
        self.assertTrue(data.shape[1] > 0)

    def test_query_activitydb(self):
        simple_query = f"SELECT * from metadata_info"
        data = query_activitydb(simple_query)
        self.assertTrue(data.shape[0] > 0)
        self.assertTrue(data.shape[1] > 0)

    def test_query_activities(self):
        self.assertTrue(True)

    def test_query_or_readin_data(self):
        self.assertTrue(True)


if __name__ == "__main__":
    unittest.main()
