short.orf <- dataset.clean %>%
  filter(relative_length_in_orthogroup < 0.9 | distance_length_in_orthogroup < -120) %>%
  select(query, new_orthogroup, mean_length_in_orthogroup, median_length_in_orthogroup, nucleotides, distance_length_in_orthogroup, relative_length_in_orthogroup) %>%
  mutate(diff = median_length_in_orthogroup - mean_length_in_orthogroup)

dataset.clean <- dataset.clean %>%
  filter(!query %in% short.orf$query) 

#delete sequences based on alignments (late start position, misalignment, missing domains etc...)
dataset.clean <- dataset.clean %>%
  filter(!(new_orthogroup == "143" & sample %in% c("AU008_S337", "CAM19_1_C025_CR_S25_L001", "CAM19_2_C108_CR_S11_L004"))) %>%
  filter(!(new_orthogroup == "126" & sample %in% c("ZU008_S79_L004", "AU053_S382"))) %>%
  filter(!(new_orthogroup == "136" & sample %in% c("YS-1fin_S58_L003", "AU028_S357"))) %>%
  filter(!(new_orthogroup == "137" & sample %in% c("CAM19_1_C057_CR_S57_L001"))) %>%
  filter(!(new_orthogroup == "149" & sample %in% c("ZU096_S8_L003"))) %>%
  filter(!(new_orthogroup == "152" & sample %in% c("ZU096_S8_L003", "ZU021_S92_L004", "ZU001_S72_L004"))) %>%
  filter(!(new_orthogroup == "155" & sample %in% c("YS-1fin_S58_L003"))) %>%
  filter(!(new_orthogroup == "162" & sample %in% c("CAM19_1_C064_CR_S64_L001"))) %>%
  filter(!(new_orthogroup == "17" & sample %in% c("CAM19_2_C134_CR_S36_L004"))) %>%
  filter(!(new_orthogroup == "180" & sample %in% c("AU074_S403"))) %>%
  filter(!(new_orthogroup == "181" & sample %in% c("AU023_S352", "BL937_S52_L003", "ZU112_S11_L003"))) %>%
  filter(!(new_orthogroup == "185" & sample %in% c("ZU009_S80_L004", "ZU015_S86_L004"))) %>%
  filter(!(new_orthogroup == "188" & sample %in% c("ZU106_S9_L003"))) %>%
  filter(!(new_orthogroup == "195" & sample %in% c("CAM19_1_C059_CR_S59_L001", "AU049_S378", "AU054_S383"))) %>%
  filter(!(new_orthogroup == "198" & nucleotides < 960)) %>%
  filter(!(new_orthogroup == "10" & nucleotides < 660)) %>%
  filter(!(new_orthogroup == "201" & sample %in% c("AU006_S335", "AU085_S414"))) %>%
  filter(!(new_orthogroup == "206" & sample %in% c("AU076_S405"))) %>%
  filter(!(new_orthogroup == "207" & sample %in% c("ZU011_S82_L004"))) %>%
  filter(!(new_orthogroup == "212" & sample %in% c("BL942_S57_L003", "CAM19_1_C019_CR_S19_L001", "CAM19_2_C003_CR_S5_L002", "CAM19_2_C149_CR_S49_L004"))) %>%
  filter(!(new_orthogroup == "221" & sample %in% c("CAM19_2_C166_CR_S65_L004"))) %>%
  filter(!(new_orthogroup == "232" & sample %in% c("YS-5_S62_L003", "YS-3-1_S60_L003"))) %>%
  filter(!(new_orthogroup == "236" & sample %in% c("CAM19_2_C020_CR_S21_L002"))) %>%
  filter(!(new_orthogroup == "238" & sample %in% c("BL940_S55_L003"))) %>%
  filter(!(new_orthogroup == "243" & sample %in% c("ZU096_S8_L003"))) %>%
  filter(!(new_orthogroup == "247" & sample %in% c("CAM19_2_C020_CR_S21_L002", "CAM19_2_C065_CR_S66_L002", "YS-1fin_S58_L003", "CAM19_2_C005_CR_S7_L002"))) %>%
  filter(!(new_orthogroup == "256" & sample %in% c("AU002_S331", "AU062_S391", "AU070_S399", "AU096_S425"))) %>%
  filter(!(new_orthogroup == "258" & sample %in% c("CAM19_2_C166_CR_S65_L004"))) %>%
  filter(!(new_orthogroup == "261" & sample %in% c("AU008_S337", "AU078_S407"))) %>%
  filter(!(new_orthogroup == "262" & nucleotides < 1050)) %>%
  filter(!(new_orthogroup == "263" & sample %in% c("YS-1fin_S58_L003"))) %>%
  filter(!(new_orthogroup == "267" & sample %in% c("CAM19_2_C166_CR_S65_L004"))) %>%
  filter(!(new_orthogroup == "28" & sample %in% c("CAM19_1_C058_CR_S58_L001", "CAM19_2_C102_CR_S6_L004", "CAM19_2_C119_CR_S22_L004", "CAM19_2_C142_CR_S43_L004", "CAM19_2_C144_CR_S45_L004", "CAM19_2_C169_CR_S68_L004", "ZU059_S1_L003"))) %>%
  filter(!(new_orthogroup == "292" & nucleotides < 1410)) %>%
  filter(!(new_orthogroup == "296" & nucleotides < 555)) %>%
  filter(!(new_orthogroup == "298" & sample %in% c("AU076_S405", "CAM19_1_C025_CR_S25_L001", "CAM19_1_C050_CR_S50_L001"))) %>%
  filter(!(new_orthogroup == "302" & sample %in% c("CAM19_1_C010_CR_S10_L001", "YS-3-1_S60_L003"))) %>%
  filter(!(new_orthogroup == "303" & nucleotides < 540)) %>%
  filter(!(new_orthogroup == "306" & sample %in% c("AU006_S335"))) %>%
  filter(!(new_orthogroup == "308" & nucleotides < 960)) %>%
  filter(!(new_orthogroup == "318" & sample %in% c("CAM19_1_C025_CR_S25_L001", "CAM19_2_C108_CR_S11_L004"))) %>%
  filter(!(new_orthogroup == "320" & nucleotides < 960)) %>%
  filter(!(new_orthogroup == "321" & nucleotides < 2040)) %>%
  filter(!(new_orthogroup == "328" & nucleotides < 1320)) %>%
  filter(!(new_orthogroup == "338" & nucleotides < 630)) %>%
  filter(!(new_orthogroup == "34" & sample %in% c("CAM19_2_C149_CR_S49_L004"))) %>%
  filter(!(new_orthogroup == "357" & sample %in% c("CAM19_2_C166_CR_S65_L004"))) %>%
  filter(!(new_orthogroup == "364" & sample %in% c("BL940_S55_L003"))) %>%
  filter(!(new_orthogroup == "372" & sample %in% c("YS-1fin_S58_L003"))) %>%
  filter(!(new_orthogroup == "376" & sample %in% c("CAM19_2_C166_CR_S65_L004"))) %>%
  filter(!(new_orthogroup == "382" & nucleotides < 1620)) %>%
  filter(!(new_orthogroup == "392" & sample %in% c("OKE010_S23_L003", "ZU106_S9_L003", "AU073_S402"))) %>%
  filter(!(new_orthogroup == "414" & sample %in% c("CAM19_1_C010_CR_S10_L001"))) %>%
  filter(!(new_orthogroup == "421" & sample %in% c("YS-3-1_S60_L003"))) %>%
  filter(!(new_orthogroup == "423" & nucleotides < 780)) %>%
  filter(!(new_orthogroup == "426" & sample %in% c("ZU008_S79_L004"))) %>%
  filter(!(new_orthogroup == "429" & sample %in% c("OKE010_S23_L003", "ZU007_S78_L004", "ZU096_S8_L003", "ZU008_S79_L004"))) %>%
  filter(!(new_orthogroup == "432" & sample %in% c("ZU021_S92_L004"))) %>%
  filter(!(new_orthogroup == "433" & sample %in% c("AU027_S356"))) %>%
  filter(!(new_orthogroup == "435" & sample %in% c("CAM19_2_C100_CR_S4_L004", "CAM19_2_C171_CR_S70_L004"))) %>%
  filter(!(new_orthogroup == "439" & sample %in% c("BL940_S55_L003"))) %>%
  filter(!(new_orthogroup == "441" & nucleotides < 480)) %>%
  filter(!(new_orthogroup == "44" & nucleotides < 660)) %>%
  filter(!(new_orthogroup == "450" & sample %in% c("ZU007_S78_L004", "AU080_S409"))) %>%
  filter(!(new_orthogroup == "451" & sample %in% c("AU076_S405"))) %>%
  filter(!(new_orthogroup == "453" & sample %in% c("AU032_S361"))) %>%
  filter(!(new_orthogroup == "480" & sample %in% c("CAM19_1_C025_CR_S25_L001"))) %>%
  filter(!(new_orthogroup == "490" & sample %in% c("YS-3-1_S60_L003"))) %>%
  filter(!(new_orthogroup == "495" & sample %in% c("ZU004_S75_L004"))) %>%
  filter(!(new_orthogroup == "499" & nucleotides < 450)) %>%
  filter(!(new_orthogroup == "503" & nucleotides < 390)) %>%
  filter(!(new_orthogroup == "504" & nucleotides < 390)) %>%
  filter(!(new_orthogroup == "520" & sample %in% c("ZU016_S87_L004", "AU080_S409"))) %>%
  filter(!(new_orthogroup == "545" & sample %in% c("CAM19_1_C050_CR_S50_L001", "AU074_S403", "YS-1fin_S58_L003", "YS-3-1_S60_L003", "YS-5_S62_L003"))) %>%
  filter(!(new_orthogroup == "568" & nucleotides < 210)) %>%
  filter(!(new_orthogroup == "578" & nucleotides < 150)) %>%
  filter(!(new_orthogroup == "57" & sample %in% c("CAM19_1_C010_CR_S10_L001"))) %>%
  filter(!(new_orthogroup == "585" & sample %in% c("AU049_S378", "AU054_S383"))) %>%
  filter(!(new_orthogroup == "58" & sample %in% c("ZU008_S79_L004", "OKE011_S24_L003"))) %>%
  filter(!(new_orthogroup == "594" & nucleotides < 390)) %>%
  filter(!(new_orthogroup == "63" & sample %in% c("CAM19_2_C166_CR_S65_L004", "OKE010_S23_L003"))) %>%
  filter(!(new_orthogroup == "68" & sample %in% c("CAM19_2_C134_CR_S36_L004"))) %>%
  filter(!(new_orthogroup == "6" & nucleotides < 1200)) %>%
  filter(!(new_orthogroup == "77" & sample %in% c("AU053_S382", "BL937_S52_L003", "CAM19_2_C133_CR_S35_L004", "ZU011_S82_L004"))) %>%
  filter(!(new_orthogroup == "78" & nucleotides < 750)) %>%
  filter(!(new_orthogroup == "81" & nucleotides < 1020)) %>%
  filter(!(new_orthogroup == "87" & nucleotides < 720)) %>%
  filter(!(new_orthogroup == "89" & nucleotides < 1860)) %>%
  filter(!(new_orthogroup == "91" & nucleotides < 546)) %>%
  filter(!(new_orthogroup == "94" & nucleotides < 1410)) %>%
  filter(!(new_orthogroup == "9" & sample %in% c("AU032_S361")))
  
  

