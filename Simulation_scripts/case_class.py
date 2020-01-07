class Case():
    def __init__(self, case_id, level):

        self.children = []

        self.case_id = case_id

        self.level = level

        #print("case ID " + str(self.case_id) + " is level " + str(self.level))

        if self.case_id != 0 and self.level == None:
            print("ERROR " + str(self.case_id))