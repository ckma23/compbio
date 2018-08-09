import json

class Energyknockdown(object):
    def energyknockdownretriever(self):
        energy_knockdown_hash = {
        "hb":
          {
          "CAT_1":.8,
          "CAT_2":.8,
          "CAT_3":.8,
          "CAT_4":.8,
          "CAT_5":.8,
          "CAT_6":.8,
          "CAT_7":.8,
          "CAT_8":.8,
          "CAT_9":.8
        },
          "vdw":
          {
          "CAT_1":.2,
          "CAT_2":.2,
          "CAT_3":.2,
          "CAT_4":.2,
          "CAT_5":.2,
          "CAT_6":.2,
          "CAT_7":.2,
          "CAT_8":.2,
          "CAT_9":.2
          }
        }

        print energy_knockdown_hash
        return energy_knockdown_hash
