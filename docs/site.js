document.body.classList.add("js-ready");

const revealTargets = document.querySelectorAll(
  ".section-card, .description-block, .metric-card, .subfigure"
);

if ("IntersectionObserver" in window) {
  const revealObserver = new IntersectionObserver(
    (entries, observer) => {
      entries.forEach((entry) => {
        if (!entry.isIntersecting) {
          return;
        }
        entry.target.classList.add("is-visible");
        observer.unobserve(entry.target);
      });
    },
    { rootMargin: "0px 0px -8% 0px", threshold: 0.08 }
  );

  revealTargets.forEach((target) => revealObserver.observe(target));
} else {
  revealTargets.forEach((target) => target.classList.add("is-visible"));
}

const navLinks = Array.from(document.querySelectorAll(".side-nav a[href^='#']"));
const sections = navLinks
  .map((link) => document.querySelector(link.getAttribute("href")))
  .filter(Boolean);

if ("IntersectionObserver" in window && sections.length) {
  const navObserver = new IntersectionObserver(
    (entries) => {
      const visible = entries
        .filter((entry) => entry.isIntersecting)
        .sort((a, b) => b.intersectionRatio - a.intersectionRatio)[0];

      if (!visible) {
        return;
      }

      const id = `#${visible.target.id}`;
      navLinks.forEach((link) => {
        link.classList.toggle("is-active", link.getAttribute("href") === id);
      });
    },
    { rootMargin: "-18% 0px -68% 0px", threshold: [0.08, 0.24, 0.5] }
  );

  sections.forEach((section) => navObserver.observe(section));
}

function closeModal(modal) {
  const image = modal.querySelector("[data-modal-img]");
  modal.hidden = true;
  if (image) {
    image.removeAttribute("src");
  }
}

document.addEventListener("click", (event) => {
  const tab = event.target.closest("[data-tab-target]");
  if (tab) {
    const tabList = tab.closest(".tab-list");
    const tabs = tabList ? tabList.querySelectorAll("[data-tab-target]") : [];
    const tabsRoot = tab.closest(".field-tabs");
    const panels = tabsRoot ? tabsRoot.querySelectorAll("[data-tab-panel]") : [];
    const target = tab.getAttribute("data-tab-target");

    tabs.forEach((button) => {
      const active = button === tab;
      button.classList.toggle("is-active", active);
      button.setAttribute("aria-selected", String(active));
    });

    panels.forEach((panel) => {
      const active = panel.getAttribute("data-tab-panel") === target;
      panel.classList.toggle("is-active", active);
      panel.hidden = !active;
    });
    return;
  }

  const openButton = event.target.closest("[data-modal-src]");
  const modal = document.querySelector("[data-modal]");
  if (openButton && modal) {
    const image = modal.querySelector("[data-modal-img]");
    const caption = modal.querySelector("[data-modal-caption]");
    image.src = openButton.getAttribute("data-modal-src");
    image.alt = openButton.getAttribute("data-modal-caption") || "";
    caption.textContent = openButton.getAttribute("data-modal-caption") || "";
    modal.hidden = false;
    return;
  }

  if (event.target.closest("[data-modal-close]") || event.target === modal) {
    closeModal(modal);
  }
});

document.addEventListener("keydown", (event) => {
  if (event.key !== "Escape") {
    return;
  }
  const modal = document.querySelector("[data-modal]");
  if (modal && !modal.hidden) {
    closeModal(modal);
  }
});
