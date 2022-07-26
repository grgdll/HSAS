function out = hsas_rd_digital_counts(fn)
# 
# out = rd_hsas(fn)
#
# read hyperOCR data from file processed with satcon to digital counts


#Instrument	InsFile	LogFile	Fitted	Immersed	Hex	CheckSum	TimeTags	F-TimeTags
#SATHED0258	HED258G.cal	2015-303-110108.raw	YES	NO	NO	NO	YES	YES

#Index	SN	INTTIME(ES)	SAMPLE(DELAY)	ES(305.66)	ES(308.98)	ES(312.29)	ES(315.61)	ES(318.93)	ES(322.24)	ES(325.56)	ES(328.88)	ES(332.21)	ES(335.53)	ES(338.85)	ES(342.17)	ES(345.50)	ES(348.82)	ES(352.14)	ES(355.47)	ES(358.80)	ES(362.12)	ES(365.45)	ES(368.78)	ES(372.11)	ES(375.44)	ES(378.77)	ES(382.10)	ES(385.43)	ES(388.76)	ES(392.09)	ES(395.42)	ES(398.75)	ES(402.09)	ES(405.42)	ES(408.76)	ES(412.09)	ES(415.43)	ES(418.76)	ES(422.10)	ES(425.43)	ES(428.77)	ES(432.11)	ES(435.44)	ES(438.78)	ES(442.12)	ES(445.46)	ES(448.80)	ES(452.14)	ES(455.48)	ES(458.81)	ES(462.15)	ES(465.50)	ES(468.84)	ES(472.18)	ES(475.52)	ES(478.86)	ES(482.20)	ES(485.54)	ES(488.88)	ES(492.22)	ES(495.57)	ES(498.91)	ES(502.25)	ES(505.59)	ES(508.94)	ES(512.28)	ES(515.62)	ES(518.96)	ES(522.31)	ES(525.65)	ES(528.99)	ES(532.34)	ES(535.68)	ES(539.02)	ES(542.37)	ES(545.71)	ES(549.05)	ES(552.39)	ES(555.74)	ES(559.08)	ES(562.42)	ES(565.77)	ES(569.11)	ES(572.45)	ES(575.79)	ES(579.14)	ES(582.48)	ES(585.82)	ES(589.16)	ES(592.50)	ES(595.84)	ES(599.19)	ES(602.53)	ES(605.87)	ES(609.21)	ES(612.55)	ES(615.89)	ES(619.23)	ES(622.57)	ES(625.91)	ES(629.25)	ES(632.58)	ES(635.92)	ES(639.26)	ES(642.60)	ES(645.94)	ES(649.27)	ES(652.61)	ES(655.95)	ES(659.28)	ES(662.62)	ES(665.95)	ES(669.29)	ES(672.62)	ES(675.95)	ES(679.29)	ES(682.62)	ES(685.95)	ES(689.28)	ES(692.61)	ES(695.94)	ES(699.28)	ES(702.60)	ES(705.93)	ES(709.26)	ES(712.59)	ES(715.92)	ES(719.24)	ES(722.57)	ES(725.89)	ES(729.22)	ES(732.54)	ES(735.87)	ES(739.19)	ES(742.51)	ES(745.83)	ES(749.15)	ES(752.47)	ES(755.79)	ES(759.11)	ES(762.43)	ES(765.75)	ES(769.06)	ES(772.38)	ES(775.69)	ES(779.00)	ES(782.32)	ES(785.63)	ES(788.94)	ES(792.25)	ES(795.56)	ES(798.87)	ES(802.18)	ES(805.48)	ES(808.79)	ES(812.09)	ES(815.40)	ES(818.70)	ES(822.00)	ES(825.30)	ES(828.60)	ES(831.90)	ES(835.20)	ES(838.50)	ES(841.79)	ES(845.09)	ES(848.38)	ES(851.68)	ES(854.97)	ES(858.26)	ES(861.55)	ES(864.84)	ES(868.12)	ES(871.41)	ES(874.70)	ES(877.98)	ES(881.26)	ES(884.54)	ES(887.82)	ES(891.10)	ES(894.38)	ES(897.66)	ES(900.93)	ES(904.21)	ES(907.48)	ES(910.75)	ES(914.02)	ES(917.29)	ES(920.56)	ES(923.83)	ES(927.09)	ES(930.36)	ES(933.62)	ES(936.88)	ES(940.14)	ES(943.40)	ES(946.65)	ES(949.91)	ES(953.16)	ES(956.42)	ES(959.67)	ES(962.92)	ES(966.17)	ES(969.41)	ES(972.66)	ES(975.90)	ES(979.14)	ES(982.39)	ES(985.63)	ES(988.86)	ES(992.10)	ES(995.33)	ES(998.57)	ES(1001.80)	ES(1005.03)	ES(1008.26)	ES(1011.48)	ES(1014.71)	ES(1017.93)	ES(1021.15)	ES(1024.37)	ES(1027.59)	ES(1030.81)	ES(1034.02)	ES(1037.24)	ES(1040.45)	ES(1043.66)	ES(1046.87)	ES(1050.07)	ES(1053.28)	ES(1056.48)	ES(1059.68)	ES(1062.88)	ES(1066.08)	ES(1069.27)	ES(1072.46)	ES(1075.66)	ES(1078.85)	ES(1082.03)	ES(1085.22)	ES(1088.40)	ES(1091.59)	ES(1094.77)	ES(1097.95)	ES(1101.12)	ES(1104.30)	ES(1107.47)	ES(1110.64)	ES(1113.81)	ES(1116.97)	ES(1120.14)	ES(1123.30)	ES(1126.46)	ES(1129.62)	ES(1132.78)	ES(1135.93)	ES(1139.08)	ES(1142.23)	DARK_SAMP(ES)	DARK_AVE(ES)	TEMP(PCB)	FRAME(COUNTER)	TIMER	CHECK(SUM)	DATETAG	TIMETAG2
#		sec	sec	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm	uW/cm^2/nm			Celsius		sec	
#1	0258	0.2560000000	0.0000000000	-0.9207583331	-0.8264237826	-0.7271381681	-0.4778115736	-0.4071159411	-0.3561771029	-0.2737434363	-0.2480580187	-0.2708502116	-0.2510742559	-0.2574123790	-0.2487287192	-0.2613095166	-0.3023014021	-0.2832371742	-0.2943309786	-0.2966543277	-0.2978119326	-0.2587658226	-0.2368259861	-0.2050015219	-0.1852701638	-0.1741640293	-0.1484036437	-0.1414150943	-0.1420635874	-0.1305061234	-0.1369269208	-0.1444862708	-0.1259376138	-0.1342971431	-0.1302681462	-0.1411262059	-0.1519744689	-0.1303749964	-0.1450856359	-0.1332298313	-0.1460752087	-0.1380559271	-0.1316912830	-0.1438354892	-0.1465777752	-0.1452615391	-0.1425778635	-0.1333642627	-0.1299466526	-0.1267302769	-0.1285818738	-0.1144820804	-0.1176512846	-0.1132454003	-0.1086166714	-0.0971945662	-0.0941149166	-0.0873338317	-0.0870514078	-0.0854941505	-0.0850023049	-0.0885697892	-0.0930705419	-0.0943776356	-0.0961808270	-0.0941151639	-0.0890987338	-0.0876669684	-0.0928465562	-0.0865206333	-0.0938209680	-0.0954204086	-0.1066921011	-0.1113923675	-0.1047634301	-0.1133118039	-0.1151242501	-0.1186672090	-0.1174431679	-0.1175568336	-0.1242955898	-0.1351557550	-0.1279000428	-0.1242815519	-0.1180030512	-0.1389711747	-0.1312938401	-0.1425415315	-0.1455766997	-0.1411285193	-0.1439839014	-0.1385663786	-0.1340937165	-0.1471738295	-0.1283657914	-0.1323159617	-0.1291781498	-0.1228140213	-0.1274261354	-0.1166903106	-0.1179661310	-0.1212152703	-0.1313842781	-0.1372490265	-0.1321426565	-0.1390927831	-0.1361986975	-0.1320267386	-0.1290184777	-0.1350233262	-0.1350568293	-0.1391583835	-0.1305633962	-0.1326681363	-0.1292783907	-0.1314030784	-0.1263298481	-0.1231006447	-0.1231339815	-0.1161167915	-0.1317425677	-0.1349227394	-0.1330814055	-0.1403566662	-0.1494364374	-0.1575431494	-0.1434492783	-0.1449591036	-0.1519627493	-0.1598417671	-0.1788927971	-0.1599684369	-0.1699083868	-0.1652267698	-0.1895754279	-0.1895850693	-0.1849102493	-0.1908148623	-0.2072711295	-0.2223050585	-0.2255615007	-0.2167307519	-0.2171753494	-0.2162517981	-0.2391503179	-0.2274323593	-0.2252641208	-0.2551814664	-0.2702899950	-0.2716851667	-0.2806827138	-0.3249467140	-0.3289246987	-0.3453784872	-0.3305378112	-0.3674892677	-0.3299992305	-0.3671788075	-0.3534237101	-0.3727253754	-0.3919249334	-0.3462296190	-0.3614044137	-0.3731510997	-0.4214059360	-0.3964957838	-0.4260649360	-0.4608194099	-0.4824452410	-0.5382661506	-0.5737231976	-0.5560752457	-0.5844072535	-0.5853052169	-0.5593567460	-0.5623201782	-0.5480248509	-0.6050533484	-0.6025415840	-0.5908339208	-0.6137130240	-0.6640694869	-0.7162193280	-0.7414398685	-0.8045134906	-0.8092412426	-0.8187739995	-0.8849176508	-0.8641620604	-0.8863602988	-0.9416546107	-1.0032975573	-0.9493077090	-0.8968656780	-1.0299506812	-1.0765390896	-1.1438705707	-1.2277461177	-1.3445135402	-1.2738767617	-1.3296174694	-1.3487048036	-1.5028594533	-1.4814772089	-1.6083752791	-1.5774188117	-1.7383914036	-1.8347904699	-1.9694114946	-1.9856965739	-2.2689374129	-2.3223272620	-2.4889166422	-2.8155393691	-2.7429981875	-3.3823078806	-3.2239228614	-3.7255737672	-3.5732546109	-4.0922860378	-4.1575128562	-4.5147326663	-5.0728356431	-5.6721227844	-6.7843877187	-7.4419500773	-7.3346555553	-7.8235124554	-8.9583257011	-8.9176735166	-10.1312319878	-11.4797138390	-12.1055324057	-12.7456444105	-12.2098347488	-12.2858979667	-13.1644847404	-14.8008748921	-15.3541844978	-16.8515265860	-16.6561320919	-17.3612355323	-19.4727500886	-22.3151661291	-22.8454480495	-22.3227600087	-22.0274603359	-21.3317907569	-17.8822473388	-16.8690858111	-14.8019783189	-11.3660434071	-10.2196112234	-9.1601205820	-8.1581995591	-7.9200037024	-7.7837659226	-8.5056792692	15	0	30.0000000000	183	7486.6800000000	124	2015-303	11:01:08.515
#2	0258	0.2560000000	0.0000000000	-0.9827213698	-0.8374256035	-0.7224937281	-0.4996354328	-0.4747578298	-0.3537073383	-0.3064071189	-0.2774362483	-0.2523980393	-0.2601111651	-0.2592162373	-0.2703954557	-0.2686602020	-0.2760479151	-0.2643169817	-0.2697843577	-0.2836452923	-0.2437776168	-0.2570517217	-0.2352242954	-0.2272105877	-0.1825543689	-0.1816068305	-0.1631919079	-0.1593121411	-0.1469950837	-0.1239365193	-0.1269773483	-0.1311716909	-0.1285909243	-0.1271904401	-0.1230674804	-0.1283283791	-0.1361642360	-0.1464552603	-0.1335833789	-0.1400158956	-0.1480279945	-0.1341469599	-0.1414300337	-0.1341675663	-0.1208086275	-0.1217471517	-0.1287349044	-0.1162172226	-0.1229073843	-0.1258758272	-0.1327134963	-0.1184758178	-0.1153403034	-0.1072971690	-0.1000051286	-0.1020513804	-0.0968035493	-0.0866809954	-0.0889639110	-0.0817447634	-0.0800723130	-0.0763390958	-0.0808567203	-0.0741513611	-0.0782689083	-0.0834990283	-0.0966977190	-0.1050642764	-0.1092446311	-0.1032397255	-0.1088409525	-0.1107308773	-0.1066921011	-0.1186076935	-0.1003600018	-0.1014215590	-0.0947796171	-0.1072109284	-0.1104570774	-0.1064757508	-0.1097268549	-0.1235305964	-0.1185265180	-0.1242815519	-0.1287267830	-0.1325947205	-0.1340665602	-0.1509119938	-0.1343805016	-0.1271065299	-0.1327684549	-0.1320311535	-0.1443349951	-0.1434669836	-0.1357324456	-0.1505853385	-0.1364167444	-0.1290729719	-0.1327417593	-0.1360727660	-0.1241156910	-0.1133029890	-0.1296233780	-0.1416620366	-0.1224431422	-0.1276694071	-0.1204476194	-0.1215855777	-0.1229613941	-0.1229629033	-0.1341970547	-0.1400168401	-0.1288491446	-0.1258286030	-0.1301304760	-0.1203682213	-0.1254829516	-0.1307241820	-0.1333362004	-0.1358763799	-0.1378498391	-0.1411611995	-0.1458737408	-0.1459920091	-0.1426742580	-0.1545671370	-0.1536276782	-0.1647656777	-0.1658073639	-0.1674595503	-0.1944513983	-0.1883957476	-0.2037049986	-0.1903597750	-0.1969740708	-0.2022841831	-0.1849102493	-0.2124231700	-0.1905760237	-0.2008279689	-0.2064603144	-0.2242727806	-0.2032596013	-0.2479252092	-0.2521372071	-0.2473888391	-0.2577229278	-0.2727723026	-0.2829794108	-0.2698142500	-0.2806827138	-0.2850874447	-0.2857811056	-0.3009888975	-0.3001159552	-0.3094471399	-0.3208503747	-0.3249686865	-0.3486161400	-0.3456184338	-0.3792837432	-0.3774093043	-0.4202838105	-0.4200529473	-0.4670181925	-0.4672875411	-0.4291174363	-0.4513286763	-0.4268010171	-0.5111915802	-0.4724222428	-0.4767420459	-0.5361805690	-0.6119860572	-0.5554459261	-0.6023415424	-0.6053170931	-0.6677501453	-0.6878898352	-0.6430308041	-0.6447614988	-0.6186782342	-0.6372617739	-0.6270980693	-0.6675008957	-0.7085392463	-0.7461255475	-0.8153465495	-0.8089027253	-0.8578127862	-0.8884504454	-0.9112282590	-1.0453513050	-0.9371419664	-1.0370209734	-1.1061660237	-1.0898743876	-1.2277461177	-1.2120352842	-1.3168425935	-1.2404538057	-1.2373819003	-1.4061920285	-1.6631594305	-1.5767612301	-1.7209891920	-1.7035562275	-1.8837639243	-1.9564396493	-2.0957204174	-2.1227571202	-1.9794125909	-2.4555815234	-2.3154879759	-2.7813894510	-3.0936394074	-3.3124879129	-3.7495058156	-3.8822572148	-4.2035703556	-4.8475813825	-5.2607228623	-6.2694998713	-6.1650277444	-6.6610979804	-7.0873350922	-7.2865638683	-7.7714466473	-8.0047215476	-8.9786323895	-10.1972236112	-10.1243336125	-11.1119161636	-11.9195607751	-13.0104679835	-13.4221835017	-15.8746722203	-16.9200992577	-15.4661444256	-17.4398192809	-18.5082265149	-19.7208293044	-18.5040942428	-21.0146643149	-18.6745047683	-19.9565569126	-18.4997741751	-21.1902318768	-18.7914128896	-16.4059706420	-15.0985358725	-11.7043613017	-10.6653925506	-9.7701175131	-9.6130780588	-10.0265637461	-8.5384343169	-7.6136369497	15	0	28.0000000000	184	7490.2300000000	0	2015-303	11:01:12.015



    fid = fopen(fn, "r");

        # skip header
        for irec = 1:4
#             disp(["----------irec = ", num2str(irec)]);
#             fflush(stdout);
            tmplab = fgets(fid);
        endfor     

        # determine if the delimiter in the extracted Satlantic file is a SPACE or a TAB
        if strcmp(tmplab(6), " ")
            DELIM = " ";
        else
            DELIM = "\t";
        endif
        

		# figure out how many channels we have
		    if (nargin == 1)
				
				tmllbl = strsplit(tmplab, DELIM);
				istart = find(cellfun(@isempty, strfind(tmllbl, "SAMPLE(DELAY)"))==0)+1;
				istop = find(cellfun(@isempty, strfind(tmllbl, "DARK_SAMP("))==0)-1;
				
		        used_pixels = 1:istop-istart+1;	
				
		    endif




        # determine if we need to read DATETAG and TIMETAG2
        
        
        
        
        ncols = length(strsplit(tmplab, DELIM));
        ncols = ncols + 3; # these are needed to account for the formated date at the end
        

        # extract wavelengths
#         out.wv = cellfun(@str2num, strsplit(tmplab, {"(", [")" DELIM], DELIM}) (8:2:end-11) )';
        out.wv = cellfun(@str2num, strsplit(tmplab, {"(", [")" DELIM], DELIM}) (8:2:end-12) )';
        out.wv = out.wv(used_pixels);
        # which instrument (ES, Li, Lt)
        instru = strsplit(tmplab, {"(", [")" DELIM], DELIM}){9:2:10};
        out.instru = instru;
        
        # save name of satcon file that was read
        out.satcon_file = fn;
        
        
        # skip another line
        tmpunits = fgets(fid);

        # now read data
        #ncols = 266;
#         fmt = [repmat(["%f" DELIM], 1, ncols), DELIM];
        fmt = [repmat(["%f" DELIM], 1, ncols-5), "%f-%f" DELIM "%f:%f:%f"];
        tmp = fscanf(fid, fmt, [ncols, Inf])';


	if isempty(tmp)
	   
  	   out.sn = [];
	   out.time = [];
	   out.int_time_sec = [];
	   out.sample_delay_sec = [];
	   out.(instru) = [];

	  disp(["no data for " fn]);
	  fflush(stdout);

	   return;

	endif



        # serial number        
        out.sn = tmp(1,2);

        # time
        out.time = datenum([tmp(:,end-4) tmp(:,end-4)*0+1 tmp(:,end-3) tmp(:,end-2) tmp(:,end-1) tmp(:,end)]);            

        # integration time (seconds)
        out.int_time_sec = tmp(:,3);
		if any(out.int_time_sec>4)
# 			disp('int time is likely in msec: converting to secs...' );
# 			fflush(stdout);
			out.int_time_sec = out.int_time_sec/1000;
		endif
			

        # sample delay time (seconds)
        out.sample_delay_sec = tmp(:,4);

        # radiometric data
        out.(instru) = tmp(:,istart:istop);   # [uW/cm^2/nm(/sr)]
        out.(instru) = out.(instru)(:,used_pixels);

#         disp("converting to Trios units");
#             fflush(stdout); 
#             out.(instru) = out.(instru)*10;   # [mW/m^2/nm(/sr)]  CONVERT TO TRIOS UNITS

     

        # DARK_SAMPLE
        out.DARK_SAMPLE = tmp(:,end-10);

        out.DARK_AVE = tmp(:,end-9);

        out.TEMP_degC = tmp(:,end-8);

        out.FRAME_count = tmp(:,end-7);

        out.TIMER = tmp(:,end-6);

        out.CHECK_SUM = tmp(:,end-5);
    


        


    fclose(fid);



endfunction

