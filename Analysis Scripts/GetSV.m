function SV = GetSV

global P

EDV = max(Get('Cavity','V','Lv'));
ESV = min(Get('Cavity','V','Lv'));

SV  = EDV - ESV;

end