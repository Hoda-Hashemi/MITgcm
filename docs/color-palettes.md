# Scientific Website Color Palettes

The live stylesheet uses **Caspian Lapis**. Each palette has light and dark values and keeps body text/status colors at WCAG AA contrast against its background or surface.

## Palette Table

| Palette | Mode | Primary | Secondary | Accent | Background | Surface/card | Text | Muted text | Success | Warning | Error | Professionalism | Scientific feel | Readability | Modern appearance |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Caspian Lapis | Light | `#155C63` | `#5E4B7B` | `#9A6A2F` | `#F4F7F8` | `#FFFFFF` | `#101820` | `#66787F` | `#1D7357` | `#8A6518` | `#9B3A44` | 10 | 10 | 10 | 9 |
| Caspian Lapis | Dark | `#72C0C4` | `#B6A2CF` | `#D0A45A` | `#081315` | `#101D20` | `#EDF6F7` | `#8EA4A9` | `#74C69D` | `#E0B65A` | `#E48D96` | 10 | 10 | 9 | 9 |
| Abyssal Copper | Light | `#17425A` | `#6C567A` | `#A1743A` | `#F6F8F6` | `#FFFFFF` | `#10212B` | `#62737A` | `#24745A` | `#89681A` | `#9A3F46` | 9 | 10 | 10 | 9 |
| Abyssal Copper | Dark | `#80B9D3` | `#C1A9CB` | `#D4A86A` | `#071117` | `#101D25` | `#EEF6F8` | `#91A4AC` | `#78CFA6` | `#E1C06A` | `#E2939B` | 9 | 10 | 9 | 9 |
| Safavid Graphite | Light | `#1E5368` | `#7A5A69` | `#8F7138` | `#F7F8F4` | `#FFFFFF` | `#141B20` | `#68737A` | `#2A7455` | `#836A1A` | `#963F4D` | 10 | 9 | 10 | 8 |
| Safavid Graphite | Dark | `#8CC7D8` | `#D1AEBE` | `#D7BA70` | `#0B1115` | `#151D22` | `#F0F4F2` | `#9AA6AB` | `#82CCA4` | `#DEC86C` | `#E69AA5` | 10 | 9 | 9 | 8 |
| Turbulence Jade | Light | `#145E56` | `#485D86` | `#98713C` | `#F3F7F5` | `#FFFFFF` | `#111C1E` | `#61767A` | `#1E7554` | `#856816` | `#963B45` | 9 | 10 | 10 | 10 |
| Turbulence Jade | Dark | `#78CABB` | `#9CB3E0` | `#D3B06D` | `#071412` | `#10211F` | `#EDF7F4` | `#8FA6A2` | `#70D19B` | `#DFC15E` | `#E48C98` | 9 | 10 | 9 | 10 |
| Persian Night Lab | Light | `#274C77` | `#6B557C` | `#9B6F35` | `#F5F6FA` | `#FFFFFF` | `#121A24` | `#65717F` | `#2B7357` | `#89681D` | `#984250` | 10 | 9 | 10 | 9 |
| Persian Night Lab | Dark | `#91B9E3` | `#C4A8D5` | `#D5AA62` | `#090F19` | `#121B29` | `#EEF3FA` | `#96A4B5` | `#7BCCA2` | `#E0BE64` | `#E4939E` | 10 | 9 | 9 | 9 |

## Preview CSS Variable Block

```css
:root {
  --bg: #f4f7f8;
  --surface: #ffffff;
  --surface-2: #edf3f4;
  --primary: #155c63;
  --secondary: #5e4b7b;
  --accent: #9a6a2f;
  --text: #101820;
  --text-soft: #344851;
  --text-muted: #66787f;
  --success: #1d7357;
  --warning: #8a6518;
  --error: #9b3a44;
}
```

## Usage

| Token | Use |
|---|---|
| Primary | Headers, active navigation, primary actions, important links. |
| Secondary | Supporting scientific/technology accents, secondary data states, subtle borders. |
| Accent | Focus rings, hover states, data links, restrained Persian-gold highlights. |
| Background | Full page background. |
| Surface/card | Cards, panels, tables, modals, and repeated experiment blocks. |
| Text | Main body copy, table data, labels that must be read quickly. |
| Muted text | Captions, units, metadata, secondary descriptions. |
| Success/warning/error | Verified, pending, and issue states. |
