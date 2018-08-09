import json

class Energyknockdown(object):
    def energyknockdownretriever(self):
        energy_knockdown_hash = {
        "hb":
          {
          "CAT_1":.9,
          "CAT_2":.9,
          "CAT_3":.9,
          "CAT_4":.9,
          "CAT_5":.9,
          "CAT_6":.9,
          "CAT_7":.9,
          "CAT_8":.9,
          "CAT_9":.9
        },
          "vdw":
          {
          "CAT_1":.1,
          "CAT_2":.1,
          "CAT_3":.1,
          "CAT_4":.1,
          "CAT_5":.1,
          "CAT_6":.1,
          "CAT_7":.1,
          "CAT_8":.1,
          "CAT_9":.1
          }
        }

        print energy_knockdown_hash
        return energy_knockdown_hash
