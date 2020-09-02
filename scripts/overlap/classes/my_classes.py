class Overlap:

    def overlap(self, interval1, interval2):
        """
        Given [0, 4] and [1, 10] returns [1, 4]
        """
        if interval2[0] <= interval1[0] <= interval2[1]:
            start = interval1[0]
        elif interval1[0] <= interval2[0] <= interval1[1]:
            start = interval2[0]
        else:
            # start = 0.0
            raise Exception("Intervals are not overlapping")

        if interval2[0] <= interval1[1] <= interval2[1]:
            end = interval1[1]
        elif interval1[0] <= interval2[1] <= interval1[1]:
            end = interval2[1]
        else:
            # end = 0.0
            raise Exception("Intervals are not overlapping")

        return (start, end)

    def percentage_overlap(self, interval1, interval2):
        """
        Given [0, 4] and [1, 10] returns 0.75
        """
        try:
            overlap = self.overlap(interval1, interval2)
        except Exception:
            return 0.0

        return (float(overlap[1]) - float(overlap[0])) / (float(interval1[1]) - float(interval1[0]))

    def reciprocal_percentage_overlap(self, interval1, interval2):
        """
        Given [0, 4] and [1, 10] returns 0.75
        """
        try:
            self.overlap(interval1, interval2)
            self.overlap(interval2, interval1)
        except Exception:
            return 0.0

        return min(self.percentage_overlap(interval1, interval2), self.percentage_overlap(interval2, interval1))