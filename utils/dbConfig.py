from neo4j import GraphDatabase

USERNAME = "neo4j"
PASSWORD = "blob1234"
driver = GraphDatabase.driver("bolt://localhost:7687", auth=(USERNAME, PASSWORD))
